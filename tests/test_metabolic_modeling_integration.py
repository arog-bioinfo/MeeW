import copy
import csv
import json
import os
import stat
import sys
import time

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

from test_annotation_integration import (
    ROOT,
    SCHEMA,
    SNAKEFILE,
    _rule_count,
    _run,
    _seed_checkpoint,
    _write_config,
)


def _write_metabolic_config(
    workdir, modules_root, *, micom=False, metabolic_qa=None, **kwargs
):
    config_path = _write_config(workdir, modules_root, **kwargs)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["metabolic_modeling"]["solver"] = "scip"
    config["metabolic_modeling"]["execution"]["mock"] = True
    config["metabolic_modeling"]["community"] = {
        "enabled": micom,
        "tools": {"smetana": False, "micom": micom, "steadycom": False},
    }
    if metabolic_qa is not None:
        config["metabolic_modeling"]["qa"]["enabled"] = metabolic_qa
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


def _seed_bakta(workdir, mag="prok_pass"):
    bakta = (
        workdir
        / "results"
        / "annotation_sample"
        / "annotation"
        / "batch"
        / "prok"
        / mag
        / "bakta"
    )
    bakta.mkdir(parents=True)
    (bakta / "annotation.faa").write_text(
        f">{mag}_protein\nMPEPTIDE\n", encoding="utf-8"
    )
    (bakta / ".snakemake_timestamp").touch()
    return bakta


def _seed_metabolic_qa(workdir, mags=("prok_pass",)):
    qa = workdir / "results" / "annotation_sample" / "metabolic_modeling" / "qa"
    qa.mkdir(parents=True)
    (qa / "filtered_mags.txt").write_text("\n".join(mags) + "\n", encoding="utf-8")
    (qa / "external_qa_summary.tsv").write_text(
        "mag\tfaa\tcompleteness\tcontamination\tstatus\n", encoding="utf-8"
    )


def _write_fake_native_carveme(root, *, fail_makedb=False):
    site = root / "fake site"
    package = site / "carveme"
    generated = package / "data" / "generated"
    generated.mkdir(parents=True)
    (package / "config.cfg").write_text(
        "[generated]\n"
        "fasta_file = data/generated/bigg_proteins.faa\n"
        "diamond_db = data/generated/bigg_proteins.dmnd\n",
        encoding="utf-8",
    )
    (generated / "bigg_proteins.faa").write_text(">bigg\nMPEPTIDE\n", encoding="utf-8")
    metadata = site / "carveme-1.6.1.dist-info"
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        "Metadata-Version: 2.1\nName: carveme\nVersion: 1.6.1\n",
        encoding="utf-8",
    )

    fake_bin = root / "fake bin"
    fake_bin.mkdir()
    diamond_calls = root / "diamond-calls.txt"
    diamond = fake_bin / "diamond"
    diamond.write_text(
        f"#!{sys.executable}\n"
        "import os, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "with Path(os.environ['MAIN_DIAMOND_CALLS']).open('a') as handle:\n"
        "    handle.write(args[0] + '\\n')\n"
        "if args[0] == 'makedb':\n"
        "    if os.environ['MAIN_DIAMOND_FAIL'] == '1':\n"
        "        raise SystemExit(29)\n"
        "    Path(args[args.index('-d') + 1] + '.dmnd').write_bytes(b'VALID-DIAMOND-INDEX\\n')\n"
        "elif args[0] == 'dbinfo':\n"
        "    database = Path(args[args.index('--db') + 1])\n"
        "    if not database.read_bytes().startswith(b'VALID'):\n"
        "        raise SystemExit(31)\n"
        "else:\n"
        "    raise SystemExit(37)\n",
        encoding="utf-8",
    )
    diamond.chmod(diamond.stat().st_mode | stat.S_IXUSR)

    carve_calls = root / "carve-calls.txt"
    carve = fake_bin / "carve"
    carve.write_text(
        f"#!{sys.executable}\n"
        "import os, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "output = Path(args[args.index('--output') + 1])\n"
        "output.parent.mkdir(parents=True, exist_ok=True)\n"
        "output.write_text('fake model\\n', encoding='utf-8')\n"
        "with Path(os.environ['MAIN_CARVE_CALLS']).open('a') as handle:\n"
        "    handle.write(output.name + '\\n')\n",
        encoding="utf-8",
    )
    carve.chmod(carve.stat().st_mode | stat.S_IXUSR)
    return site, fake_bin, diamond_calls, carve_calls, fail_makedb


def _configure_fake_native_environment(monkeypatch, fixture):
    site, fake_bin, diamond_calls, carve_calls, fail_makedb = fixture
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv(
        "PYTHONPATH",
        f"{site}{os.pathsep}{os.environ.get('PYTHONPATH', '')}",
    )
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")
    monkeypatch.setenv("MAIN_DIAMOND_CALLS", str(diamond_calls))
    monkeypatch.setenv("MAIN_CARVE_CALLS", str(carve_calls))
    monkeypatch.setenv(
        "MAIN_DIAMOND_FAIL",
        "1" if fail_makedb else "0",
    )


def _metabolic_target(kind, mag="prok_pass"):
    root = "results/annotation_sample/metabolic_modeling"
    targets = {
        "model": f"{root}/gems/{mag}.xml",
        "report": f"{root}/memote/{mag}/report.html",
        "score": f"{root}/memote/{mag}/score.json",
        "memote_summary": f"{root}/memote/memote_summary.tsv",
        "memote_quality": f"{root}/memote/memote_quality.tsv",
        "filtered_models": f"{root}/memote/filtered_models.txt",
        "qa_summary": f"{root}/qa/external_qa_summary.tsv",
        "qa_filtered": f"{root}/qa/filtered_mags.txt",
        "medium_sum": f"{root}/media/summed_demand.csv",
        "medium_m9": f"{root}/media/M9.csv",
        "micom": f"{root}/micom/summed_demand/fluxes.tsv",
        "micom_m9": f"{root}/micom/M9/fluxes.tsv",
    }
    return targets[kind]


def test_real_annotation_and_metabolic_modules_share_one_namespace_safely(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "namespace"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)

    listed = _run(workdir, config_path, "--list-rules")
    assert listed.returncode == 0, listed.stdout + listed.stderr
    rules = set(listed.stdout.splitlines())
    assert "annotation_qa_selection" in rules
    assert "metabolic_modeling_external_qa_filter" in rules
    assert "metabolic_modeling_carveme" in rules
    assert "metabolic_modeling_memote" in rules
    assert "metabolic_modeling_all" not in rules
    assert "metabolic_modeling_stage_checkm2_faa" not in rules
    assert "metabolic_modeling_checkm2" not in rules

    rendered = listed.stdout + listed.stderr
    assert "metabolic_modeling.done" not in rendered

    source = SNAKEFILE.read_text(encoding="utf-8")
    assert source.count("use rule * from metabolic_modeling") == 1
    assert source.count("metabolic_modeling.bind_protein_sources(") == 1
    for removed_plumbing in (
        "ANNOTATION_IMPORTED_NAMESPACE",
        "METABOLIC_MODELING_NAMESPACE",
        "METABOLIC_MODELING_MEMOTE_NAMESPACE",
        "METABOLIC_MODELING_IMPORTED_NAMESPACE",
        "METABOLIC_MODELING_QA_CHECKPOINT",
        "METABOLIC_MODELING_QA_BINDINGS",
        "METABOLIC_MODELING_MEMOTE_CHECKPOINT",
        "FilterCheckpointBindings",
        "bind_filter_checkpoints",
        "use rule external_qa_filter from metabolic_modeling",
        "use rule filter_memote from metabolic_modeling",
        "use rule carveme from metabolic_modeling",
    ):
        assert removed_plumbing not in source

    harness = workdir / "inspect-native-checkpoints.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("QA_CHECKPOINT_NATIVE=" + str(\n'
        "    metabolic_modeling.checkpoints.external_qa_filter is not None\n"
        "))\n"
        'print("MEMOTE_CHECKPOINT_NATIVE=" + str(\n'
        "    metabolic_modeling.checkpoints.filter_memote is not None\n"
        "))\n",
        encoding="utf-8",
    )
    inspected = _run(workdir, config_path, "--list-rules", snakefile=harness)
    assert inspected.returncode == 0, inspected.stdout + inspected.stderr
    assert "QA_CHECKPOINT_NATIVE=True" in inspected.stdout
    assert "MEMOTE_CHECKPOINT_NATIVE=True" in inspected.stdout


def test_direct_bakta_output_object_wires_both_qa_and_carveme(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "direct-dag"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)
    _seed_checkpoint(workdir, [("prok_pass", "prok", "pass")])
    bakta = _seed_bakta(workdir)

    dry_run = _run(
        workdir,
        config_path,
        "--dry-run",
        "--printshellcmds",
        _metabolic_target("report"),
    )
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    rendered = dry_run.stdout + dry_run.stderr
    assert str(bakta.relative_to(workdir)) in rendered
    assert "annotation.faa" not in rendered
    qa_job = rendered.split("checkpoint metabolic_modeling_external_qa_filter:", 1)[
        1
    ].split("rule metabolic_modeling_carveme:", 1)[0]
    carveme_job = rendered.split("rule metabolic_modeling_carveme:", 1)[1].split(
        "rule metabolic_modeling_memote:", 1
    )[0]
    assert str(bakta.relative_to(workdir)) in qa_job
    assert "annotation/qa/prok/checkm2/quality_report.tsv" in qa_job
    assert str(bakta.relative_to(workdir)) in carveme_job
    assert "<TBD>" in carveme_job
    assert "gtdbtk/classification.tsv" not in carveme_job
    assert _rule_count(rendered, "metabolic_modeling_initialize_carveme_index") == 0


def test_native_main_partial_resume_plans_only_missing_or_science_stale_models(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "native-partial-resume"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["metabolic_modeling"]["execution"]["mock"] = False
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    mags = ("cached", "missing")
    _seed_checkpoint(workdir, [(mag, "prok", "pass") for mag in mags])
    bakta = {mag: _seed_bakta(workdir, mag) for mag in mags}
    _seed_metabolic_qa(workdir, mags)
    cached = workdir / _metabolic_target("model", "cached")
    cached.parent.mkdir(parents=True)
    cached.write_text("cached model\n", encoding="utf-8")

    missing = _run(
        workdir,
        config_path,
        "--dry-run",
        "--printshellcmds",
        "--rerun-triggers=mtime",
        *(_metabolic_target("model", mag) for mag in mags),
    )
    assert missing.returncode == 0, missing.stdout + missing.stderr
    rendered = missing.stdout + missing.stderr
    assert _rule_count(rendered, "metabolic_modeling_initialize_carveme_index") == 1
    assert _rule_count(rendered, "metabolic_modeling_carveme") == 1
    assert rendered.index("rule metabolic_modeling_initialize_carveme_index:") < (
        rendered.index("rule metabolic_modeling_carveme:")
    )
    assert _metabolic_target("model", "missing") in rendered
    assert f"output: {_metabolic_target('model', 'cached')}" not in rendered

    newer = cached.stat().st_mtime_ns + 1_000_000_000
    os.utime(bakta["cached"] / ".snakemake_timestamp", ns=(newer, newer))
    science_stale = _run(
        workdir,
        config_path,
        "--dry-run",
        "--printshellcmds",
        "--rerun-triggers=mtime",
        *(_metabolic_target("model", mag) for mag in mags),
    )
    assert science_stale.returncode == 0, science_stale.stdout + science_stale.stderr
    rendered = science_stale.stdout + science_stale.stderr
    assert _rule_count(rendered, "metabolic_modeling_initialize_carveme_index") == 1
    assert _rule_count(rendered, "metabolic_modeling_carveme") == 2
    assert f"output: {_metabolic_target('model', 'cached')}" in rendered


def test_native_main_partial_resume_executes_one_preparation_and_one_missing_model(
    tmp_path, all_modules_root, monkeypatch
):
    workdir = tmp_path / "native-main-success"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["metabolic_modeling"]["execution"]["mock"] = False
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    mags = ("cached", "missing")
    _seed_checkpoint(workdir, [(mag, "prok", "pass") for mag in mags])
    for mag in mags:
        _seed_bakta(workdir, mag)
    _seed_metabolic_qa(workdir, mags)
    cached = workdir / _metabolic_target("model", "cached")
    cached.parent.mkdir(parents=True)
    cached.write_text("preserved model\n", encoding="utf-8")
    cached_before = (cached.read_bytes(), cached.stat().st_mtime_ns)
    fixture = _write_fake_native_carveme(workdir / "native fixture")
    _configure_fake_native_environment(monkeypatch, fixture)

    result = _run(
        workdir,
        config_path,
        "--rerun-triggers=mtime",
        *(_metabolic_target("model", mag) for mag in mags),
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (cached.read_bytes(), cached.stat().st_mtime_ns) == cached_before
    assert (workdir / _metabolic_target("model", "missing")).read_text() == (
        "fake model\n"
    )
    _, _, diamond_calls, carve_calls, _ = fixture
    assert diamond_calls.read_text().splitlines() == ["makedb", "dbinfo"]
    assert carve_calls.read_text().splitlines() == ["missing.xml"]
    rendered = result.stdout + result.stderr
    assert _rule_count(rendered, "metabolic_modeling_initialize_carveme_index") == 1
    assert _rule_count(rendered, "metabolic_modeling_carveme") == 1


def test_mock_partial_resume_keeps_cached_score_and_builds_missing_cohort_member(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "partial-memote-resume"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)
    mags = ("cached", "missing")
    _seed_checkpoint(workdir, [(mag, "prok", "pass") for mag in mags])
    for mag in mags:
        _seed_bakta(workdir, mag)
    _seed_metabolic_qa(workdir, mags)

    cached_model = workdir / _metabolic_target("model", "cached")
    cached_model.parent.mkdir(parents=True)
    cached_model.write_text("cached model\n", encoding="utf-8")
    cached_report = workdir / _metabolic_target("report", "cached")
    cached_score = workdir / _metabolic_target("score", "cached")
    cached_report.parent.mkdir(parents=True)
    cached_report.write_text("<html>cached report</html>\n", encoding="utf-8")
    cached_score.write_text(
        json.dumps(
            {
                "score": {
                    "total_score": 1.0,
                    "sections": [{"section": "consistency", "score": 1.0}],
                }
            }
        ),
        encoding="utf-8",
    )
    cached_qa_summary = workdir / _metabolic_target("qa_summary")
    cached_qa_filtered = workdir / _metabolic_target("qa_filtered")
    protected = {
        path: (path.read_bytes(), path.stat().st_mtime_ns)
        for path in (
            cached_model,
            cached_report,
            cached_score,
            cached_qa_summary,
            cached_qa_filtered,
        )
    }

    harness = workdir / "partial-metabolic-targets.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        "rule integrated_metabolic_targets:\n"
        "    input:\n"
        "        expected_metabolic_modeling_outputs,\n",
        encoding="utf-8",
    )
    result = _run(
        workdir,
        config_path,
        "--rerun-triggers=mtime",
        "--printshellcmds",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in protected
    } == protected
    assert (workdir / _metabolic_target("score", "missing")).is_file()

    rendered = result.stdout + result.stderr
    assert "localcheckpoint metabolic_modeling_external_qa_filter:" not in rendered
    assert rendered.count("localrule metabolic_modeling_carveme:") == 1
    assert rendered.count("localrule metabolic_modeling_memote:") == 1
    with (workdir / _metabolic_target("memote_quality")).open(
        newline="", encoding="utf-8"
    ) as handle:
        quality_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["mag"], row["status"]) for row in quality_rows] == [
        ("cached", "PASS"),
        ("missing", "PASS"),
    ]


def test_failed_native_main_preparation_blocks_missing_model_and_preserves_cached(
    tmp_path, all_modules_root, monkeypatch
):
    workdir = tmp_path / "native-main-failure"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["metabolic_modeling"]["execution"]["mock"] = False
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    mags = ("cached", "missing")
    _seed_checkpoint(workdir, [(mag, "prok", "pass") for mag in mags])
    for mag in mags:
        _seed_bakta(workdir, mag)
    _seed_metabolic_qa(workdir, mags)
    cached = workdir / _metabolic_target("model", "cached")
    cached.parent.mkdir(parents=True)
    cached.write_text("preserved model\n", encoding="utf-8")
    cached_before = (cached.read_bytes(), cached.stat().st_mtime_ns)
    fixture = _write_fake_native_carveme(
        workdir / "native failure fixture", fail_makedb=True
    )
    _configure_fake_native_environment(monkeypatch, fixture)

    result = _run(
        workdir,
        config_path,
        "--rerun-triggers=mtime",
        *(_metabolic_target("model", mag) for mag in mags),
    )
    assert result.returncode != 0
    assert (cached.read_bytes(), cached.stat().st_mtime_ns) == cached_before
    assert not (workdir / _metabolic_target("model", "missing")).exists()
    _, _, diamond_calls, carve_calls, _ = fixture
    assert diamond_calls.read_text().splitlines() == ["makedb"]
    assert not carve_calls.exists()
    assert not list(
        (workdir / ".snakemake" / "meew-metabolic-modeling" / "carveme-index").glob(
            "*.ready.json"
        )
    )
    rendered = result.stdout + result.stderr
    assert _rule_count(rendered, "metabolic_modeling_initialize_carveme_index") == 1
    assert "localrule metabolic_modeling_carveme:" not in rendered


def test_enabled_gtdbtk_classification_is_a_direct_carveme_dependency(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "gtdb-dependency"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root, gtdbtk=True)
    _seed_checkpoint(workdir, [("prok_pass", "prok", "pass")])
    _seed_bakta(workdir)
    metabolic_qa = (
        workdir / "results" / "annotation_sample" / "metabolic_modeling" / "qa"
    )
    metabolic_qa.mkdir(parents=True)
    (metabolic_qa / "filtered_mags.txt").write_text("prok_pass\n", encoding="utf-8")
    (metabolic_qa / "external_qa_summary.tsv").write_text(
        "mag\tfaa\tcompleteness\tcontamination\tstatus\n",
        encoding="utf-8",
    )

    dry_run = _run(
        workdir,
        config_path,
        "--dry-run",
        "--printshellcmds",
        _metabolic_target("model"),
    )
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    rendered = dry_run.stdout + dry_run.stderr
    assert _rule_count(rendered, "annotation_gtdbtk") == 1
    assert _rule_count(rendered, "metabolic_modeling_carveme") == 1
    carveme_job = rendered.split("rule metabolic_modeling_carveme:", 1)[1]
    assert (
        "results/annotation_sample/annotation/gtdbtk/classification.tsv" in carveme_job
    )
    assert "Input files updated by another job" in carveme_job


def test_enabled_gtdbtk_mixed_domains_drive_mock_carveme_without_source_changes(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "gtdb-mixed"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root, gtdbtk=True)
    mags = ["EST6_prok_bin_002", "EST6_prok_bin_005"]
    _seed_checkpoint(
        workdir,
        [(mag, "prok", "pass") for mag in mags],
    )
    bakta_dirs = [_seed_bakta(workdir, mag) for mag in mags]
    classification = (
        workdir
        / "results"
        / "annotation_sample"
        / "annotation"
        / "gtdbtk"
        / "classification.tsv"
    )
    classification.parent.mkdir(parents=True)
    classification.write_text(
        "domain\tuser_genome\tclassification\n"
        "bacteria\tEST6_prok_bin_002\td__Bacteria\n"
        "archaea\tEST6_prok_bin_005\td__Archaea\n",
        encoding="utf-8",
    )
    protected = [
        classification,
        *(bakta / "annotation.faa" for bakta in bakta_dirs),
    ]
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in protected}

    result = _run(
        workdir,
        config_path,
        *(_metabolic_target("model", mag) for mag in mags),
    )
    assert result.returncode == 0, result.stdout + result.stderr

    for mag, expected_universe in zip(mags, ("bacteria", "archaea")):
        log = (
            workdir
            / "results"
            / "annotation_sample"
            / "metabolic_modeling"
            / "logs"
            / "carveme"
            / f"{mag}.log"
        ).read_text(encoding="utf-8")
        assert f"MAG={mag} universe={expected_universe} source=GTDB" in log
        assert f"source_path={classification.relative_to(workdir)}" in log
    assert {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in protected
    } == before


@pytest.mark.parametrize(
    ("rows", "expected_jobs"),
    [
        ([("prok_pass", "prok", "pass")], True),
        ([("euk_pass", "euk", "pass")], False),
        ([("prok_fail", "prok", "fail")], False),
        ([], False),
    ],
    ids=["qualifying-prok", "euk-only", "all-filtered", "empty"],
)
def test_checkpoint_aware_main_targets_include_only_qualifying_prok_routes(
    tmp_path, all_modules_root, rows, expected_jobs
):
    workdir = tmp_path / "route"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)
    _seed_checkpoint(workdir, rows)
    if expected_jobs:
        _seed_bakta(workdir)
        qa_result = _run(
            workdir,
            config_path,
            _metabolic_target("qa_summary"),
            _metabolic_target("qa_filtered"),
        )
        assert qa_result.returncode == 0, qa_result.stdout + qa_result.stderr

    harness = workdir / "metabolic-targets.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        "rule integrated_metabolic_targets:\n"
        "    input:\n"
        "        expected_metabolic_modeling_outputs,\n",
        encoding="utf-8",
    )
    dry_run = _run(
        workdir,
        config_path,
        "--dry-run",
        "--printshellcmds",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    rendered = dry_run.stdout + dry_run.stderr
    assert (
        _rule_count(rendered, "metabolic_modeling_filter_memote") > 0
    ) is expected_jobs
    # Once first QA is cached, its finite cohort is resolved before the
    # downstream Memote checkpoint and all score producers are concrete.
    for rule in (
        "metabolic_modeling_carveme",
        "metabolic_modeling_memote",
    ):
        assert (_rule_count(rendered, rule) == 1) is expected_jobs
    assert _rule_count(rendered, "metabolic_modeling_memote_summary") == 0
    assert _rule_count(rendered, "metabolic_modeling_external_qa_filter") == 0
    if expected_jobs:
        assert "metabolic_modeling/memote/memote_quality.tsv" in rendered
        assert "metabolic_modeling/memote/filtered_models.txt" in rendered
        assert _metabolic_target("model") in rendered
        assert _metabolic_target("report") in rendered
        assert _metabolic_target("score") in rendered
        assert "metabolic_modeling/memote/memote_summary.tsv" not in rendered
        filter_job = rendered.split("checkpoint metabolic_modeling_filter_memote:", 1)[
            1
        ].split("\n\n", 1)[0]
        assert _metabolic_target("score") in filter_job
        assert _metabolic_target("qa_summary") not in filter_job
        assert _metabolic_target("qa_filtered") not in filter_job
        assert "<TBD>" not in filter_job
    else:
        assert "/metabolic_modeling/" not in rendered
    assert "metabolic_modeling.done" not in rendered


@pytest.mark.parametrize(
    ("annotation_qa", "requested_metabolic_qa"),
    [(True, True), (True, False), (False, True), (False, False)],
    ids=[
        "annotation-on-metabolic-on",
        "annotation-on-metabolic-off",
        "annotation-off-metabolic-on",
        "annotation-off-metabolic-off",
    ],
)
def test_annotation_and_requested_metabolic_qa_policy_uses_owner_defaults(
    tmp_path,
    all_modules_root,
    annotation_qa,
    requested_metabolic_qa,
):
    workdir = tmp_path / "qa-policy"
    workdir.mkdir()
    config_path = _write_metabolic_config(
        workdir,
        all_modules_root,
        qa=annotation_qa,
        qa_filter=annotation_qa,
        metabolic_qa=requested_metabolic_qa,
    )
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    main_resources = workdir / "main-resources"
    module_resources = workdir / "module-resources" if requested_metabolic_qa else None
    config["resources_root"] = str(main_resources)
    config["metabolic_modeling"]["resources_root"] = (
        str(module_resources) if module_resources else ""
    )
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    _seed_checkpoint(
        workdir,
        [("prok_pass", "prok", "pass")],
        qa=annotation_qa,
    )
    bakta = _seed_bakta(workdir)
    effective_qa = annotation_qa and requested_metabolic_qa

    harness = workdir / "inspect-metabolic-config.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("EFFECTIVE_QA=" + str(metabolic_modeling.qa_enabled))\n'
        'print("QA_TSV=" + metabolic_modeling.qa_tsv)\n'
        'print("QA_COLUMNS=" + ",".join((\n'
        "    metabolic_modeling.qa_id_col,\n"
        "    metabolic_modeling.qa_completeness_col,\n"
        "    metabolic_modeling.qa_contamination_col,\n"
        ")))\n"
        'print("CHECKM2_ENABLED=" + str(metabolic_modeling.checkm2_enabled))\n'
        'print("QA_USES_CHECKM2=" + str(metabolic_modeling.qa_uses_checkm2))\n'
        'print("INPUT_CONFIG=" + str(metabolic_modeling.input_config))\n'
        'print("PROJECT=" + METABOLIC_MODELING_MODULE_CONFIG["sample"])\n'
        'print("COVERM=" + METABOLIC_MODELING_MODULE_CONFIG["micom"]["abundance_tsv"])\n'
        'print("RESOURCES=" + METABOLIC_MODELING_MODULE_CONFIG["resources_root"])\n'
        'print("LOGS_ROOT=" + metabolic_modeling.LOGS_ROOT)\n'
        'print("CARVEME_LOG=" + str(\n'
        '    workflow.get_rule("metabolic_modeling_carveme").log[0]\n'
        "))\n",
        encoding="utf-8",
    )

    listed = _run(workdir, config_path, "--list-rules", snakefile=harness)
    assert listed.returncode == 0, listed.stdout + listed.stderr
    rules = set(listed.stdout.splitlines())
    assert ("annotation_qa_checkm2_prok" in rules) is annotation_qa
    assert "metabolic_modeling_checkm2" not in rules
    assert "metabolic_modeling_stage_checkm2_faa" not in rules
    assert f"EFFECTIVE_QA={effective_qa}" in listed.stdout
    expected_report = (
        "results/annotation_sample/annotation/qa/prok/checkm2/quality_report.tsv"
        if effective_qa
        else ""
    )
    assert f"QA_TSV={expected_report}" in listed.stdout
    assert "QA_COLUMNS=Name,Completeness,Contamination" in listed.stdout
    assert "CHECKM2_ENABLED=False" in listed.stdout
    assert "QA_USES_CHECKM2=False" in listed.stdout
    assert "INPUT_CONFIG={}" in listed.stdout
    assert "PROJECT=annotation_sample" in listed.stdout
    assert (
        "COVERM=results/annotation_sample/binning/coverm/"
        "prok_genome_coverage.tsv" in listed.stdout
    )
    expected_resources = module_resources or main_resources
    assert f"RESOURCES={expected_resources}" in listed.stdout
    assert "LOGS_ROOT=<results>/logs" in listed.stdout
    assert (
        "CARVEME_LOG=results/{sample}/metabolic_modeling/logs/"
        "carveme/{mag}.log" in listed.stdout
    )

    qa_execution = _run(
        workdir,
        config_path,
        _metabolic_target("qa_summary"),
        _metabolic_target("qa_filtered"),
    )
    assert qa_execution.returncode == 0, qa_execution.stdout + qa_execution.stderr
    with (workdir / _metabolic_target("qa_summary")).open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {
            "mag": "prok_pass",
            "faa": str((bakta / "annotation.faa").relative_to(workdir)),
            "completeness": "90.0" if effective_qa else "NA",
            "contamination": "1.0" if effective_qa else "NA",
            "status": "PASS",
        }
    ]
    assert (
        workdir
        / "results"
        / "annotation_sample"
        / "metabolic_modeling"
        / "logs"
        / "qa"
        / "external_qa_filter.log"
    ).is_file()

    direct = _run(workdir, config_path, "--dry-run", _metabolic_target("model"))
    assert direct.returncode == 0, direct.stdout + direct.stderr
    assert _rule_count(direct.stdout + direct.stderr, "metabolic_modeling_carveme") == 1


@pytest.mark.parametrize("genome_context", ["prok", "euk"])
def test_singular_routes_do_not_evaluate_or_import_metabolic_modeling(
    tmp_path,
    all_modules_root,
    genome_context,
):
    workdir = tmp_path / genome_context
    workdir.mkdir()
    config_path = _write_metabolic_config(
        workdir,
        all_modules_root,
        genome_context=genome_context,
        include_binning=False,
    )

    listed = _run(workdir, config_path, "--list-rules")
    assert listed.returncode == 0, listed.stdout + listed.stderr
    assert not any(
        rule.startswith("metabolic_modeling_") for rule in listed.stdout.splitlines()
    )


def test_integrated_metabolic_schema_rejects_main_owned_adapter_settings():
    schema = yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema)
    shipped = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    validator.validate(shipped)

    assert shipped["metabolic_modeling"]["solver"] == "cplex"
    assert shipped["metabolic_modeling"]["resources_root"] == ""
    assert shipped["metabolic_modeling"]["simulation_medium"] == {
        "source": "discover",
        "aggregation": "sum",
        "predefined_medium": "M9",
        "min_growth": 0.1,
        "open_exchanges": True,
        "exchange_bound": 1000,
    }
    configured_root = copy.deepcopy(shipped)
    configured_root["resources_root"] = "/tmp/shared resources"
    configured_root["metabolic_modeling"]["resources_root"] = "/tmp/model resources"
    validator.validate(configured_root)
    legacy_roots = copy.deepcopy(shipped)
    legacy_roots.pop("resources_root")
    legacy_roots["metabolic_modeling"].pop("resources_root")
    validator.validate(legacy_roots)
    for key, value in (
        ("sample", "other"),
        ("input", {"faa_samples": "handoff.tsv"}),
        ("checkm2", {"enabled": True}),
        ("pathvars", {"results": "elsewhere"}),
        ("stage_sheets", {"enabled": True}),
        ("micom", {"abundance_tsv": "user-owned.tsv"}),
    ):
        invalid = copy.deepcopy(shipped)
        invalid["metabolic_modeling"][key] = value
        with pytest.raises(ValidationError):
            validator.validate(invalid)

    for key, value in (
        ("tsv", "quality.tsv"),
        (
            "columns",
            {
                "id": "bin",
                "completeness": "complete",
                "contamination": "contam",
            },
        ),
    ):
        invalid = copy.deepcopy(shipped)
        invalid["metabolic_modeling"]["qa"][key] = value
        with pytest.raises(ValidationError):
            validator.validate(invalid)


@pytest.mark.parametrize(
    ("community_enabled", "micom_enabled"),
    [(False, True), (True, False)],
    ids=["community-disabled", "micom-disabled"],
)
def test_disabled_community_or_micom_does_not_require_prok_relative_abundance(
    tmp_path, all_modules_root, community_enabled, micom_enabled
):
    workdir = tmp_path / "micom-disabled"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root, micom=False)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["metabolic_modeling"]["community"]["enabled"] = community_enabled
    config["metabolic_modeling"]["community"]["tools"]["micom"] = micom_enabled
    config["binning"]["coverm"]["domains"] = ["euk"]
    config["binning"]["coverm"]["methods"] = ["mean"]
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    Draft202012Validator(yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))).validate(
        config
    )
    listed = _run(workdir, config_path, "--list-rules")
    assert listed.returncode == 0, listed.stdout + listed.stderr


@pytest.mark.parametrize(
    ("domains", "methods", "message"),
    [
        (["euk"], ["relative_abundance"], "requires 'prok'"),
        (["prok"], ["mean"], "requires 'relative_abundance'"),
    ],
)
def test_enabled_micom_rejects_missing_prok_relative_abundance(
    tmp_path, all_modules_root, domains, methods, message
):
    workdir = tmp_path / "invalid-micom-abundance"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root, micom=True)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["binning"]["coverm"]["domains"] = domains
    config["binning"]["coverm"]["methods"] = methods
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    validator = Draft202012Validator(yaml.safe_load(SCHEMA.read_text(encoding="utf-8")))
    with pytest.raises(ValidationError):
        validator.validate(config)

    listed = _run(workdir, config_path, "--list-rules")
    assert listed.returncode != 0
    assert message in listed.stdout + listed.stderr


def test_micom_material_input_is_exact_prok_coverm_producer(tmp_path, all_modules_root):
    workdir = tmp_path / "micom-coverm-dag"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root, micom=True)
    _seed_checkpoint(workdir, [("prok_pass", "prok", "pass")])
    _seed_bakta(workdir)

    metabolic_root = workdir / "results" / "annotation_sample" / "metabolic_modeling"
    qa_root = metabolic_root / "qa"
    qa_root.mkdir(parents=True)
    (qa_root / "filtered_mags.txt").write_text("prok_pass\n", encoding="utf-8")
    (qa_root / "external_qa_summary.tsv").write_text(
        "mag\tfaa\tcompleteness\tcontamination\tstatus\n", encoding="utf-8"
    )
    memote_root = metabolic_root / "memote"
    memote_root.mkdir(parents=True)
    (memote_root / "filtered_models.txt").write_text("prok_pass\n", encoding="utf-8")
    (memote_root / "memote_quality.tsv").write_text(
        "mag\ttotal_score\tconsistency\tannotation\tstatus\n", encoding="utf-8"
    )

    result = _run(workdir, config_path, "--d3dag", _metabolic_target("micom"))
    assert result.returncode == 0, result.stdout + result.stderr
    graph = json.loads(result.stdout)
    coverm_ids = {
        node["id"]
        for node in graph["nodes"]
        if node["value"]["rule"] == "binning_coverm_genome_coverage"
    }
    micom_ids = {
        node["id"]
        for node in graph["nodes"]
        if node["value"]["rule"] == "metabolic_modeling_micom"
    }
    assert coverm_ids and micom_ids
    assert any(
        edge["u"] in coverm_ids and edge["v"] in micom_ids for edge in graph["edges"]
    )

    dry_run = _run(
        workdir,
        config_path,
        "--dry-run",
        "--printshellcmds",
        _metabolic_target("micom"),
    )
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    micom_job = (
        (dry_run.stdout + dry_run.stderr)
        .split("rule metabolic_modeling_micom:", 1)[1]
        .split("\n\n[", 1)[0]
    )
    assert (
        "results/annotation_sample/binning/coverm/prok_genome_coverage.tsv" in micom_job
    )
    assert "binning/coverm/euk_genome_coverage.tsv" not in micom_job


def test_predefined_m9_schedules_one_shared_medium_and_zero_discovery_jobs(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "predefined-m9"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root, micom=True)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["metabolic_modeling"]["simulation_medium"]["source"] = "predefined"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    _seed_checkpoint(workdir, [("prok_pass", "prok", "pass")])
    _seed_bakta(workdir)
    _seed_metabolic_qa(workdir)
    metabolic_root = workdir / "results" / "annotation_sample" / "metabolic_modeling"
    model = metabolic_root / "gems" / "prok_pass.xml"
    model.parent.mkdir(parents=True)
    model.write_text("cached model\n", encoding="utf-8")
    memote = metabolic_root / "memote"
    memote.mkdir(parents=True)
    (memote / "filtered_models.txt").write_text("prok_pass\n", encoding="utf-8")
    (memote / "memote_quality.tsv").write_text(
        "mag\ttotal_score\tconsistency\tannotation\tstatus\n",
        encoding="utf-8",
    )

    dry_run = _run(
        workdir,
        config_path,
        "--dry-run",
        "--printshellcmds",
        _metabolic_target("micom_m9"),
    )
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    rendered = dry_run.stdout + dry_run.stderr
    assert _rule_count(rendered, "metabolic_modeling_prepare_simulation_medium") == 1
    assert _rule_count(rendered, "metabolic_modeling_discover_individual_medium") == 0
    assert _rule_count(rendered, "metabolic_modeling_micom") == 1
    assert _metabolic_target("medium_m9") in rendered


def test_mock_direct_route_executes_unconditional_no_community_final_targets(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "mock-execution"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)
    _seed_checkpoint(workdir, [("prok_pass", "prok", "pass")])
    bakta = _seed_bakta(workdir)
    targets = [
        _metabolic_target("qa_summary"),
        _metabolic_target("qa_filtered"),
        _metabolic_target("model"),
        _metabolic_target("report"),
        _metabolic_target("score"),
        _metabolic_target("memote_summary"),
        _metabolic_target("memote_quality"),
        _metabolic_target("filtered_models"),
    ]
    harness = workdir / "metabolic-targets.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        "rule integrated_metabolic_targets:\n"
        "    input:\n"
        "        expected_metabolic_modeling_outputs,\n",
        encoding="utf-8",
    )

    result = _run(
        workdir,
        config_path,
        "--printshellcmds",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for target in targets:
        assert (workdir / target).is_file()

    with (workdir / _metabolic_target("qa_summary")).open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {
            "mag": "prok_pass",
            "faa": str((bakta / "annotation.faa").relative_to(workdir)),
            "completeness": "90.0",
            "contamination": "1.0",
            "status": "PASS",
        }
    ]
    assert (
        json.loads((workdir / _metabolic_target("score")).read_text(encoding="utf-8"))[
            "score"
        ]["total_score"]
        == 1.0
    )
    with (workdir / _metabolic_target("memote_summary")).open(
        newline="", encoding="utf-8"
    ) as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["mag"], row["total_score"]) for row in summary_rows] == [
        ("prok_pass", "1.0")
    ]
    with (workdir / _metabolic_target("memote_quality")).open(
        newline="", encoding="utf-8"
    ) as handle:
        quality_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["mag"], row["status"]) for row in quality_rows] == [
        ("prok_pass", "PASS")
    ]
    assert (workdir / _metabolic_target("filtered_models")).read_text(
        encoding="utf-8"
    ) == "prok_pass\n"

    cached = _run(
        workdir,
        config_path,
        "--printshellcmds",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert cached.returncode == 0, cached.stdout + cached.stderr
    assert "Nothing to be done" in cached.stdout + cached.stderr


def test_main_medium_default_cache_and_older_custom_replacement_use_mtime_only(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "medium-default-cache"
    workdir.mkdir()
    config_path = _write_metabolic_config(workdir, all_modules_root)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["metabolic_modeling"]["community"] = {
        "enabled": True,
        "tools": {"smetana": False, "micom": False, "steadycom": True},
    }
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    _seed_checkpoint(workdir, [("prok_pass", "prok", "pass")])
    _seed_bakta(workdir)

    harness = workdir / "metabolic-medium-targets.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        "rule integrated_metabolic_targets:\n"
        "    input:\n"
        "        expected_metabolic_modeling_outputs,\n",
        encoding="utf-8",
    )

    fresh = _run(
        workdir,
        config_path,
        "--rerun-triggers=mtime",
        "--printshellcmds",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert fresh.returncode == 0, fresh.stdout + fresh.stderr
    rendered = fresh.stdout + fresh.stderr
    assert (
        rendered.count("localrule metabolic_modeling_discover_individual_medium:") == 1
    )
    assert (
        rendered.count("localrule metabolic_modeling_prepare_simulation_medium:") == 1
    )
    assert rendered.count("localrule metabolic_modeling_steadycom:") == 1

    cached = _run(
        workdir,
        config_path,
        "--rerun-triggers=mtime",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert cached.returncode == 0, cached.stdout + cached.stderr
    assert "Nothing to be done" in cached.stdout + cached.stderr

    custom_a = workdir / "custom-a.csv"
    custom_b = workdir / "custom-b.csv"
    custom_a.write_text("nutrient_id,uptake\nmock_nutrient_e,1\n", encoding="utf-8")
    custom_b.write_text("nutrient_id,uptake\nmock_nutrient_e,2\n", encoding="utf-8")
    os.utime(custom_b, ns=(1_000_000_000, 1_000_000_000))
    time.sleep(0.02)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    medium = config["metabolic_modeling"]["simulation_medium"]
    medium["source"] = "predefined"
    medium["predefined_medium"] = str(custom_a)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    first_custom = _run(
        workdir,
        config_path,
        "--rerun-triggers=mtime",
        "--printshellcmds",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert first_custom.returncode == 0, first_custom.stdout + first_custom.stderr
    rendered = first_custom.stdout + first_custom.stderr
    assert "localrule metabolic_modeling_discover_individual_medium:" not in rendered
    assert (
        rendered.count("localrule metabolic_modeling_prepare_simulation_medium:") == 1
    )
    custom_output = (
        workdir
        / "results"
        / "annotation_sample"
        / "metabolic_modeling"
        / "media"
        / "custom.csv"
    )
    assert custom_output.read_text(encoding="utf-8") == (
        "nutrient_id,uptake\nmock_nutrient_e,1.0\n"
    )

    custom_cached = _run(
        workdir,
        config_path,
        "--rerun-triggers=mtime",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert custom_cached.returncode == 0, custom_cached.stdout + custom_cached.stderr
    assert "Nothing to be done" in custom_cached.stdout + custom_cached.stderr
    assert custom_b.stat().st_mtime_ns < custom_output.stat().st_mtime_ns

    time.sleep(0.02)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["metabolic_modeling"]["simulation_medium"]["predefined_medium"] = str(
        custom_b
    )
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    replacement = _run(
        workdir,
        config_path,
        "--rerun-triggers=mtime",
        "--printshellcmds",
        "integrated_metabolic_targets",
        snakefile=harness,
    )
    assert replacement.returncode == 0, replacement.stdout + replacement.stderr
    rendered = replacement.stdout + replacement.stderr
    assert "localrule metabolic_modeling_discover_individual_medium:" not in rendered
    assert (
        rendered.count("localrule metabolic_modeling_prepare_simulation_medium:") == 1
    )
    assert custom_output.read_text(encoding="utf-8") == (
        "nutrient_id,uptake\nmock_nutrient_e,2.0\n"
    )
