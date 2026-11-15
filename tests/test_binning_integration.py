import copy
import json
import subprocess
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"

ROUTES = {
    "short": {
        "sample": "binning_short",
        "header": ["group", "sample", "data_type", "read1", "read2"],
        "rows": [
            [
                "binning_short",
                "short",
                "short_pe",
                "inputs/short_R1.fastq",
                "inputs/short_R2.fastq",
            ]
        ],
        "assembly_rule": "assembly_assemble_spades",
        "mapping_rule": "assembly_bwamem2_map",
        "assembly_dir": "assembly-spades",
        "mapping_dir": "mapping-bwamem2",
        "semibin2_type": "short_read",
    },
    "long": {
        "sample": "binning_long",
        "header": ["group", "sample", "data_type", "read1"],
        "rows": [["binning_long", "long", "nano-hq", "inputs/long.fastq"]],
        "assembly_rule": "assembly_assemble_flye_long",
        "mapping_rule": "assembly_minimap2_map",
        "assembly_dir": "assembly-flye",
        "mapping_dir": "mapping-minimap2",
        "semibin2_type": "long_read",
    },
    "hybrid": {
        "sample": "binning_hybrid",
        "header": ["group", "sample", "data_type", "read1", "read2"],
        "rows": [
            [
                "binning_hybrid",
                "short",
                "short_pe",
                "inputs/short_R1.fastq",
                "inputs/short_R2.fastq",
            ],
            ["binning_hybrid", "long", "pacbio-hifi", "inputs/long.fastq", ""],
        ],
        "assembly_rule": "assembly_assemble_spades",
        "mapping_rule": "assembly_bwamem2_map",
        "assembly_dir": "assembly-spades",
        "mapping_dir": "mapping-bwamem2",
        "semibin2_type": "long_read",
    },
}

DIRECT_RULES = {
    "comebin": {"binning_comebin"},
    "metabinner": {
        "binning_metabinner_depth",
        "binning_metabinner_profiles",
        "binning_metabinner",
        "binning_metabinner_bins",
    },
    "semibin2": {"binning_semibin2"},
    "metabat2": {"binning_metabat2_depth", "binning_metabat2"},
}

MINIMAL_BINNER_OPTIONS = {
    "comebin": {"extra": ""},
    "metabinner": {},
    "semibin2": {"environment": "global", "extra": ""},
    "metabat2": {"min_contig_length": 1500, "extra": ""},
}


def _write_run_sheet(workdir, route):
    run_sheet = workdir / "config" / "run.tsv"
    run_sheet.parent.mkdir(parents=True)
    run_sheet.write_text(
        "\t".join(route["header"])
        + "\n"
        + "".join("\t".join(row) + "\n" for row in route["rows"]),
        encoding="utf-8",
    )
    path_columns = {
        index
        for index, column in enumerate(route["header"])
        if column in {"read1", "read2", "contam"}
    }
    for row in route["rows"]:
        for index in path_columns:
            if row[index]:
                path = workdir / row[index]
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
    return run_sheet


def _write_config(
    workdir,
    modules_root,
    route,
    tools,
    *,
    refinement=False,
    coverm_domains=None,
):
    config = yaml.safe_load(
        (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    config["modules_root"] = str(modules_root)
    config["preprocessing"]["run_sheet"] = str(_write_run_sheet(workdir, route))
    config["metabolic_modeling"]["community"]["enabled"] = False
    config["binning"]["binning"]["tools"] = tools
    config["binning"]["refinement"].update(
        {
            "enabled": refinement,
            "method": "das_tool" if refinement else "none",
            "euk_source_tool": tools[0],
        }
    )
    if coverm_domains is not None:
        config["binning"]["coverm"]["domains"] = coverm_domains
    config_path = workdir / "config" / "main.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


def _minimal_binning_fragment(tool):
    return {
        "binning": {
            "tools": [tool],
            tool: copy.deepcopy(MINIMAL_BINNER_OPTIONS[tool]),
        },
        "refinement": {"enabled": False, "method": "none"},
        "domain_separation": {
            "fallback_domain": "unknown",
            "tiara": {"extra": ""},
        },
        "coverm": {
            "input_mode": "bam",
            "domains": ["prok"],
            "mapper": "minimap2-sr",
            "methods": ["relative_abundance"],
            "threads": 1,
            "extra": "",
        },
    }


def _run(workdir, config_path, *arguments, snakefile=SNAKEFILE):
    return subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
            "--replace-workflow-config",
            "--nolock",
            "--cores",
            "1",
            *arguments,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def _domain_target(route):
    return f"results/{route['sample']}/binning/domain/domain_assignment.tsv"


def _assembly_paths(route):
    root = f"results/{route['sample']}/assembly"
    return (
        f"{root}/{route['assembly_dir']}/assembly.fasta",
        f"{root}/{route['mapping_dir']}/alignment.bam",
    )


def _graph(result):
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(result.stdout)


def _assert_summary(workdir, config_path, target):
    result = _run(workdir, config_path, "--summary", target)
    assert result.returncode == 0, result.stdout + result.stderr
    assert target in result.stdout


def _rule_ids(graph, rule):
    return {node["id"] for node in graph["nodes"] if node["value"]["rule"] == rule}


def _has_rule_edge(graph, upstream, downstream):
    upstream_ids = _rule_ids(graph, upstream)
    downstream_ids = _rule_ids(graph, downstream)
    return any(
        edge["u"] in upstream_ids and edge["v"] in downstream_ids
        for edge in graph["edges"]
    )


@pytest.mark.parametrize("route_name", ROUTES)
def test_binning_uses_exact_route_selected_assembly_materials(
    tmp_path, route_name, all_modules_root
):
    route = ROUTES[route_name]
    workdir = tmp_path / route_name
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        route,
        ["comebin"],
    )
    target = _domain_target(route)

    graph = _graph(_run(workdir, config_path, "--d3dag", target))
    assert _has_rule_edge(graph, route["assembly_rule"], "binning_comebin")
    assert _has_rule_edge(graph, route["mapping_rule"], "binning_comebin")
    assert _has_rule_edge(graph, route["assembly_rule"], "binning_tiara_predict")

    dry_run = _run(workdir, config_path, "--dry-run", "--printshellcmds", target)
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    output = dry_run.stdout + dry_run.stderr
    contigs, bam = _assembly_paths(route)
    comebin_job = output.split("rule binning_comebin:", 1)[1].split("\n\nrule ", 1)[0]
    assert contigs in comebin_job
    assert bam in comebin_job
    assert f"rule binning_tiara_predict:\n    input: {contigs}" in output
    opposite_mapping = (
        "mapping-minimap2"
        if route["mapping_dir"] == "mapping-bwamem2"
        else "mapping-bwamem2"
    )
    assert opposite_mapping not in output
    assert "binning_all" not in {node["value"]["rule"] for node in graph["nodes"]}


@pytest.mark.parametrize("tool", DIRECT_RULES)
def test_each_direct_binner_uses_its_imported_internal_dag(
    tmp_path, tool, all_modules_root
):
    route = ROUTES["short"]
    workdir = tmp_path / tool
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        route,
        [tool],
    )

    graph = _graph(_run(workdir, config_path, "--d3dag", _domain_target(route)))
    rules = {node["value"]["rule"] for node in graph["nodes"]}
    assert DIRECT_RULES[tool] <= rules
    assert rules.isdisjoint(
        set().union(
            *(ruleset for name, ruleset in DIRECT_RULES.items() if name != tool)
        )
    )
    assert {"binning_tiara_predict", "binning_separate_bins_by_domain"} <= rules

    assembly_rule = route["assembly_rule"]
    mapping_rule = route["mapping_rule"]
    if tool == "comebin":
        assert _has_rule_edge(graph, assembly_rule, "binning_comebin")
        assert _has_rule_edge(graph, mapping_rule, "binning_comebin")
    elif tool == "metabinner":
        assert _has_rule_edge(graph, mapping_rule, "binning_metabinner_depth")
        assert _has_rule_edge(graph, assembly_rule, "binning_metabinner_profiles")
        assert _has_rule_edge(
            graph, "binning_metabinner_depth", "binning_metabinner_profiles"
        )
        assert _has_rule_edge(
            graph, "binning_metabinner_profiles", "binning_metabinner"
        )
        assert _has_rule_edge(graph, "binning_metabinner", "binning_metabinner_bins")
    elif tool == "semibin2":
        assert _has_rule_edge(graph, assembly_rule, "binning_semibin2")
        assert _has_rule_edge(graph, mapping_rule, "binning_semibin2")
    else:
        assert _has_rule_edge(graph, mapping_rule, "binning_metabat2_depth")
        assert _has_rule_edge(graph, assembly_rule, "binning_metabat2")
        assert _has_rule_edge(graph, "binning_metabat2_depth", "binning_metabat2")
    _assert_summary(workdir, config_path, _domain_target(route))


def test_multi_binner_das_tool_refinement_keeps_imported_dependencies(
    tmp_path, all_modules_root
):
    route = ROUTES["short"]
    workdir = tmp_path / "refinement"
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        route,
        ["comebin", "semibin2", "metabat2", "metabinner"],
        refinement=True,
    )

    graph = _graph(_run(workdir, config_path, "--d3dag", _domain_target(route)))
    rules = {node["value"]["rule"] for node in graph["nodes"]}
    assert {
        "binning_comebin",
        "binning_semibin2",
        "binning_metabat2_depth",
        "binning_metabat2",
        "binning_metabinner_depth",
        "binning_metabinner_profiles",
        "binning_metabinner",
        "binning_metabinner_bins",
        "binning_tiara_predict",
        "binning_prepare_refinement_inputs",
        "binning_refine_das_tool",
        "binning_separate_bins_by_domain",
    } <= rules
    for upstream in (
        "binning_comebin",
        "binning_semibin2",
        "binning_metabat2",
        "binning_metabinner_bins",
    ):
        assert _has_rule_edge(graph, upstream, "binning_prepare_refinement_inputs")
    assert _has_rule_edge(
        graph, "binning_tiara_predict", "binning_prepare_refinement_inputs"
    )
    assert _has_rule_edge(
        graph, route["assembly_rule"], "binning_prepare_refinement_inputs"
    )
    assert _has_rule_edge(graph, route["assembly_rule"], "binning_refine_das_tool")
    assert _has_rule_edge(
        graph, "binning_prepare_refinement_inputs", "binning_refine_das_tool"
    )
    assert _has_rule_edge(
        graph, "binning_refine_das_tool", "binning_separate_bins_by_domain"
    )
    assert _has_rule_edge(graph, route["mapping_rule"], "binning_metabinner_depth")
    assert _has_rule_edge(graph, route["assembly_rule"], "binning_metabinner_profiles")

    _assert_summary(workdir, config_path, _domain_target(route))


@pytest.mark.parametrize("route_name", ROUTES)
def test_semibin2_type_is_derived_from_main_assembly_route(
    tmp_path, route_name, all_modules_root
):
    route = ROUTES[route_name]
    workdir = tmp_path / f"semibin-{route_name}"
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        route,
        ["semibin2"],
    )
    harness = tmp_path / f"inspect-semibin-{route_name}.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("MAIN_SEMIBIN2_TYPE=" + '
        'BINNING_MODULE_CONFIG["binning"]["semibin2"]["sequencing_type"])\n',
        encoding="utf-8",
    )

    result = _run(workdir, config_path, "--list-rules", snakefile=harness)
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"MAIN_SEMIBIN2_TYPE={route['semibin2_type']}" in result.stdout


def test_coverm_uses_final_configured_domains_and_exact_assembly_bam(
    tmp_path, all_modules_root
):
    route = ROUTES["long"]
    workdir = tmp_path / "coverm"
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        route,
        ["metabat2"],
        coverm_domains=["euk"],
    )
    target = f"results/{route['sample']}/binning/coverm/euk_genome_coverage.tsv"

    graph = _graph(_run(workdir, config_path, "--d3dag", target))
    assert len(_rule_ids(graph, "binning_coverm_genome_coverage")) == 1
    assert _has_rule_edge(
        graph, "binning_separate_bins_by_domain", "binning_coverm_genome_coverage"
    )
    assert _has_rule_edge(
        graph, route["mapping_rule"], "binning_coverm_genome_coverage"
    )

    dry_run = _run(workdir, config_path, "--dry-run", "--printshellcmds", target)
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    output = dry_run.stdout + dry_run.stderr
    _, bam = _assembly_paths(route)
    assert (
        "rule binning_coverm_genome_coverage:\n"
        f"    input: results/{route['sample']}/binning/domain/euk_bins, {bam}"
    ) in output
    assert "prok_genome_coverage.tsv" not in output
    assert "unknown_genome_coverage.tsv" not in output
    _assert_summary(workdir, config_path, target)

    all_targets = _run(workdir, config_path, "--dry-run")
    assert all_targets.returncode == 0, all_targets.stdout + all_targets.stderr
    rendered_targets = all_targets.stdout + all_targets.stderr
    assert target in rendered_targets
    assert "binning/coverm/prok_genome_coverage.tsv" not in rendered_targets


def test_final_targets_include_unknown_domain_contract(tmp_path, all_modules_root):
    route = ROUTES["short"]
    workdir = tmp_path / "targets"
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        route,
        ["comebin"],
    )

    result = _run(workdir, config_path, "--dry-run")
    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    domain_root = f"results/{route['sample']}/binning/domain"
    for target in (
        "prok_bins",
        "euk_bins",
        "unknown_bins",
        "domain_assignment.tsv",
        "annotation_samples.tsv",
        "bin_rename_map.tsv",
    ):
        assert f"{domain_root}/{target}" in output
    for domain in ("prok", "euk"):
        assert (
            f"results/{route['sample']}/binning/coverm/" f"{domain}_genome_coverage.tsv"
        ) in output
    listed = _run(workdir, config_path, "--list-rules")
    assert listed.returncode == 0, listed.stdout + listed.stderr
    listed_rules = set(listed.stdout.splitlines())
    assert "binning_all" not in listed_rules


def test_main_binning_schema_accepts_current_modes_and_rejects_owned_inputs():
    schema = yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema)
    shipped = yaml.safe_load(
        (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    validator.validate(shipped)
    assert shipped["resources_root"] == ""
    assert shipped["binning"]["resources_root"] == ""
    assert shipped["binning"]["checkm1"] == {
        "database_path": "",
        "download_database": True,
    }

    fragment = yaml.safe_load(
        (ROOT / "config" / "config_binning.yaml").read_text(encoding="utf-8")
    )
    Draft202012Validator(schema["properties"]["binning"]).validate(fragment)

    legacy = copy.deepcopy(shipped)
    legacy.pop("resources_root")
    legacy["binning"].pop("resources_root")
    legacy["binning"].pop("checkm1")
    validator.validate(legacy)

    for tool in DIRECT_RULES:
        direct = copy.deepcopy(shipped)
        direct["binning"]["binning"]["tools"] = [tool]
        validator.validate(direct)

    refinement = copy.deepcopy(shipped)
    refinement["binning"]["binning"]["tools"] = ["comebin", "metabat2"]
    refinement["binning"]["refinement"].update(
        {"enabled": True, "method": "das_tool", "euk_source_tool": "comebin"}
    )
    validator.validate(refinement)

    invalid_mutations = [
        lambda value: value["binning"].__setitem__("sample_sheet", "input.tsv"),
        lambda value: value["binning"].__setitem__("stage_sheets", {"enabled": False}),
        lambda value: value["binning"]["binning"]["semibin2"].__setitem__(
            "sequencing_type", "short_read"
        ),
        lambda value: value["binning"]["coverm"].__setitem__("input_mode", "coupled"),
        lambda value: value["binning"]["coverm"].__setitem__("bam", "user-owned.bam"),
        lambda value: value["binning"]["coverm"].__setitem__("domains", ["unknown"]),
        lambda value: value["binning"]["coverm"].__setitem__("enabled", False),
        lambda value: value["binning"]["binning"].__setitem__(
            "tools", ["comebin", "metabat2"]
        ),
    ]
    for mutate in invalid_mutations:
        invalid = copy.deepcopy(shipped)
        mutate(invalid)
        with pytest.raises(ValidationError):
            validator.validate(invalid)


@pytest.mark.parametrize("tool", MINIMAL_BINNER_OPTIONS)
def test_main_binning_schema_accepts_minimal_selected_tool_fragment(tool):
    schema = yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema["properties"]["binning"])

    validator.validate(_minimal_binning_fragment(tool))


@pytest.mark.parametrize("tool", MINIMAL_BINNER_OPTIONS)
def test_main_binning_schema_rejects_selected_tool_without_options_block(tool):
    schema = yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema["properties"]["binning"])
    fragment = _minimal_binning_fragment(tool)
    del fragment["binning"][tool]

    with pytest.raises(ValidationError):
        validator.validate(fragment)


def test_main_rejects_refinement_source_that_is_not_selected(
    tmp_path, all_modules_root
):
    route = ROUTES["short"]
    workdir = tmp_path / "invalid-euk-source"
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        route,
        ["comebin", "metabat2"],
        refinement=True,
    )
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["binning"]["refinement"]["euk_source_tool"] = "semibin2"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    result = _run(workdir, config_path, "--list-rules")

    assert result.returncode != 0
    assert "refinement.euk_source_tool must be one of the selected" in (
        result.stdout + result.stderr
    )


@pytest.mark.parametrize(
    ("main_root", "module_root", "explicit", "ambient", "expected"),
    [
        ("/tmp/main resources", "", "", "", "/tmp/main resources/checkm1_db"),
        ("/tmp/main resources", "/tmp/bin root", "", "", "/tmp/bin root/checkm1_db"),
        (
            "/tmp/main resources",
            "/tmp/bin root",
            "/tmp/exact CheckM1",
            "",
            "/tmp/exact CheckM1",
        ),
        ("", "", "", "/tmp/ambient resources", "/tmp/ambient resources/checkm1_db"),
    ],
    ids=["main", "module", "explicit", "ambient"],
)
@pytest.mark.parametrize("tool", ["comebin", "metabinner"])
def test_main_resource_root_and_checkm1_precedence(
    tmp_path,
    monkeypatch,
    all_modules_root,
    main_root,
    module_root,
    explicit,
    ambient,
    expected,
    tool,
):
    route = ROUTES["short"]
    workdir = tmp_path / "resources"
    workdir.mkdir()
    config_path = _write_config(workdir, all_modules_root, route, [tool])
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["resources_root"] = main_root
    config["binning"]["resources_root"] = module_root
    config["binning"]["checkm1"] = {
        "database_path": explicit,
        "download_database": True,
    }
    config["annotation"]["resources_root"] = module_root
    config["metabolic_modeling"]["resources_root"] = module_root
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    if ambient:
        monkeypatch.setenv("MEEW_RESOURCES", ambient)
    else:
        monkeypatch.delenv("MEEW_RESOURCES", raising=False)

    harness = workdir / "inspect-resources.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("CHECKM1_DATABASE=" + binning.CHECKM1_DATABASE_PATH)\n'
        'print("BINNING_ROOT=" + BINNING_MODULE_CONFIG.get("resources_root", ""))\n'
        'print("ANNOTATION_ROOT=" + ANNOTATION_MODULE_CONFIG.get("resources_root", ""))\n'
        'print("METABOLIC_ROOT=" + METABOLIC_MODELING_MODULE_CONFIG.get("resources_root", ""))\n',
        encoding="utf-8",
    )
    result = _run(workdir, config_path, "--list-rules", snakefile=harness)
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"CHECKM1_DATABASE={expected}" in result.stdout
    inherited = module_root or main_root
    for name in ("BINNING_ROOT", "ANNOTATION_ROOT", "METABOLIC_ROOT"):
        assert f"{name}={inherited}" in result.stdout


def test_main_legacy_blank_checkm1_defaults_to_enabled(tmp_path, all_modules_root):
    route = ROUTES["short"]
    workdir = tmp_path / "legacy-checkm1"
    workdir.mkdir()
    config_path = _write_config(workdir, all_modules_root, route, ["comebin"])
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config.pop("resources_root", None)
    config["binning"].pop("resources_root", None)
    config["binning"].pop("checkm1", None)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    harness = workdir / "inspect-legacy-checkm1.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("CHECKM1_DOWNLOAD=" + str(binning.CHECKM1_DOWNLOAD_DATABASE))\n',
        encoding="utf-8",
    )
    result = _run(workdir, config_path, "--list-rules", snakefile=harness)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "CHECKM1_DOWNLOAD=True" in result.stdout
