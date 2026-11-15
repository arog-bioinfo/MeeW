import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_PATH = ROOT / ".github" / "workflows" / "main.yaml"
DISABLED_WORKFLOW_PATH = ROOT / ".github" / "workflows-disabled" / "main.yaml"

MODULE_CHECKOUTS = {
    "arog-bioinfo/MeeW-Preprocessing": ".ci/modules/MeeW-Preprocessing",
    "arog-bioinfo/MeeW-Assembly": ".ci/modules/MeeW-Assembly",
    "arog-bioinfo/MeeW-Binning": ".ci/modules/MeeW-Binning",
    "arog-bioinfo/MeeW-Annotation": ".ci/modules/MeeW-Annotation",
    "arog-bioinfo/MeeW-Metabolic_Modeling": (".ci/modules/MeeW-Metabolic_Modeling"),
}

ROUTES = [
    {
        "name": "metagenome",
        "id": "metagenome",
        "config": ".test/config/metagenome.yaml",
    },
    {
        "name": "prokaryotic isolate",
        "id": "prokaryotic-isolate",
        "config": ".test/config/prokaryotic-isolate.yaml",
    },
    {
        "name": "eukaryotic isolate",
        "id": "eukaryotic-isolate",
        "config": ".test/config/eukaryotic-isolate.yaml",
    },
]


def _workflow():
    return yaml.safe_load(WORKFLOW_PATH.read_text(encoding="utf-8"))


def _triggers(workflow):
    # PyYAML's YAML 1.1 resolver treats the unquoted Actions key `on` as true.
    return workflow.get("on", workflow.get(True))


def _checkout_steps(job):
    return [step for step in job["steps"] if step.get("uses") == "actions/checkout@v6"]


def test_active_workflow_has_template_topology_and_controls():
    assert WORKFLOW_PATH.is_file()
    assert not DISABLED_WORKFLOW_PATH.exists()

    workflow = _workflow()
    triggers = _triggers(workflow)
    assert triggers["push"]["branches"] == ["main"]
    assert triggers["pull_request"]["branches"] == ["main"]
    assert "workflow_dispatch" in triggers
    assert workflow["concurrency"] == {
        "group": "${{ github.workflow }}-${{ github.ref }}",
        "cancel-in-progress": True,
    }
    assert workflow["permissions"] == {"contents": "read"}
    assert workflow["env"] == {
        "MEEW_MODULES_ROOT": "${{ github.workspace }}/.ci/modules"
    }

    jobs = workflow["jobs"]
    assert set(jobs) == {"formatting", "linting", "python-tests", "testing"}
    assert jobs["testing"]["needs"] == [
        "formatting",
        "linting",
        "python-tests",
    ]
    assert {name: job["timeout-minutes"] for name, job in jobs.items()} == {
        "formatting": 15,
        "linting": 15,
        "python-tests": 20,
        "testing": 20,
    }
    assert jobs["testing"]["strategy"]["fail-fast"] is False
    assert jobs["testing"]["strategy"]["matrix"]["route"] == ROUTES


def test_actions_python_and_snakemake_versions_match_module_ci():
    workflow = _workflow()
    jobs = workflow["jobs"]
    formatting = jobs["formatting"]

    assert formatting["steps"][0]["uses"] == "actions/checkout@v6"
    formatter = formatting["steps"][1]
    assert formatter["uses"] == "super-linter/super-linter@v8"
    assert formatter["env"] == {
        "DEFAULT_BRANCH": "main",
        "GITHUB_TOKEN": "${{ secrets.GITHUB_TOKEN }}",
        "VALIDATE_ALL_CODEBASE": True,
        "VALIDATE_PYTHON_BLACK": True,
        "VALIDATE_SHELL_SHFMT": True,
        "VALIDATE_SNAKEMAKE_SNAKEFMT": True,
        "VALIDATE_YAML_PRETTIER": True,
    }

    python_steps = jobs["python-tests"]["steps"]
    setup_python = next(
        step for step in python_steps if step.get("uses") == "actions/setup-python@v6"
    )
    assert setup_python["with"]["python-version"] == "3.13"
    assert any(
        step.get("run") == "python -m pip install -r tests/requirements-workflow.txt"
        for step in python_steps
    )
    assert any(
        step.get("run") == "python -m pytest -q -ra tests" for step in python_steps
    )

    snakemake_steps = [
        step
        for job in (jobs["linting"], jobs["testing"])
        for step in job["steps"]
        if step.get("uses") == "snakemake/snakemake-github-action@v2"
    ]
    assert len(snakemake_steps) == 2
    assert all(
        step["with"]["snakemake-version"] == "9.19.0" for step in snakemake_steps
    )


def test_integrated_jobs_checkout_all_public_main_modules_without_custom_tokens():
    workflow = _workflow()
    for job_name in ("linting", "python-tests", "testing"):
        checkout_steps = _checkout_steps(workflow["jobs"][job_name])
        current = [
            step for step in checkout_steps if "repository" not in step.get("with", {})
        ]
        assert len(current) == 1

        module_steps = [
            step for step in checkout_steps if "repository" in step.get("with", {})
        ]
        assert {
            step["with"]["repository"]: step["with"]["path"] for step in module_steps
        } == MODULE_CHECKOUTS
        assert all(step["with"]["ref"] == "main" for step in module_steps)
        assert all("token" not in step.get("with", {}) for step in checkout_steps)

    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    assert set(re.findall(r"secrets\.([A-Za-z0-9_]+)", source)) == {"GITHUB_TOKEN"}
    assert source.count("${{ secrets.GITHUB_TOKEN }}") == 1


def test_requirements_are_exactly_pinned():
    assert (ROOT / "tests" / "requirements.txt").read_text(
        encoding="utf-8"
    ).splitlines() == [
        "pytest==9.0.2",
        "PyYAML==6.0.3",
        "jsonschema==4.26.0",
        "pandas==2.3.3",
    ]
    assert (ROOT / "tests" / "requirements-workflow.txt").read_text(
        encoding="utf-8"
    ).splitlines() == [
        "-r requirements.txt",
        "snakemake==9.19.0",
    ]


def test_route_fixtures_are_minimal_deterministic_dag_inputs():
    expected_contexts = {
        "metagenome": "metagenome",
        "prokaryotic-isolate": "prok",
        "eukaryotic-isolate": "euk",
    }
    for route, context in expected_contexts.items():
        config_path = ROOT / ".test" / "config" / f"{route}.yaml"
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        assert config["genome_context"] == context
        assert config["preprocessing"]["auto_preprocess"] == {"enabled": False}
        assert config["annotation"]["qa"]["enabled"] is False
        assert config["annotation"]["qa_filter"]["enabled"] is False
        assert config["metabolic_modeling"]["community"]["enabled"] is False
        run_sheet = ROOT / ".test" / config["preprocessing"]["run_sheet"]
        assert run_sheet.is_file()

    assert (ROOT / ".test" / "data" / "short_R1.fastq").is_file()
    assert (ROOT / ".test" / "data" / "short_R2.fastq").is_file()
    assert not (ROOT / ".test" / "results").exists()


def test_matrix_execution_is_dry_run_only_and_has_no_heavy_ci_features():
    workflow = _workflow()
    testing = workflow["jobs"]["testing"]
    execution_steps = [
        step
        for step in testing["steps"]
        if step.get("uses") == "snakemake/snakemake-github-action@v2"
    ]
    assert len(execution_steps) == 1
    execution = execution_steps[0]
    assert execution["with"]["directory"] == ".test"
    assert execution["with"]["snakefile"] == "workflow/Snakefile"
    args = execution["with"]["args"]
    for required in ("--dry-run", "--printshellcmds", "--cores 2"):
        assert args.count(required) == 1
    assert "--configfile ${{ matrix.route.config }}" in args

    rendered_testing = yaml.safe_dump(testing).lower()
    for forbidden in (
        "--report",
        "software-deployment",
        "use-conda",
        "conda-cleanup",
        "upload-artifact",
        "download-artifact",
    ):
        assert forbidden not in rendered_testing
    assert all("run" not in step for step in testing["steps"])
