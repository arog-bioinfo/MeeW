import copy
import csv
import json
import re
import subprocess
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"
SELECTION_COLUMNS = [
    "sample",
    "path",
    "domain",
    "completeness",
    "contamination",
    "min_completeness",
    "max_contamination",
    "decision",
    "reason",
]


def _write_config(
    workdir,
    modules_root,
    *,
    genome_context="metagenome",
    gtdbtk=False,
    qa=True,
    qa_filter=True,
    include_binning=True,
):
    inputs = workdir / "inputs"
    inputs.mkdir(parents=True)
    for name in ("R1.fastq", "R2.fastq"):
        (inputs / name).touch()
    run_sheet = workdir / "run.tsv"
    run_sheet.write_text(
        "group\tsample\tdata_type\tread1\tread2\n"
        "annotation_sample\tshort\tshort_pe\tinputs/R1.fastq\tinputs/R2.fastq\n",
        encoding="utf-8",
    )
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config["modules_root"] = str(modules_root)
    config["genome_context"] = genome_context
    config["preprocessing"]["run_sheet"] = str(run_sheet)
    if not include_binning:
        config.pop("binning")
    config["annotation"]["qa"]["enabled"] = qa
    config["annotation"]["qa_filter"]["enabled"] = qa_filter
    config["annotation"]["gtdbtk"]["enabled"] = gtdbtk
    for section, field, value in (
        ("bakta", "db", workdir / "resources" / "bakta"),
        ("gtdbtk", "data_dir", workdir / "resources" / "gtdbtk"),
        ("metaeuk", "cache_dir", workdir / "resources" / "metaeuk"),
        (
            "recognizer_prok",
            "resources_dir",
            workdir / "resources" / "recognizer",
        ),
        (
            "recognizer_euk",
            "resources_dir",
            workdir / "resources" / "recognizer",
        ),
        ("upimapi", "resources_dir", workdir / "resources" / "upimapi"),
    ):
        config["annotation"][section][field] = str(value)
    if genome_context == "euk":
        config["annotation"]["funannotate2"]["species"] = "Integration fixture"
    config_path = workdir / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


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


def _graph(result):
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(result.stdout)


def _rule_ids(graph, name):
    return {node["id"] for node in graph["nodes"] if node["value"]["rule"] == name}


def _has_edge(graph, upstream, downstream):
    upstream_ids = _rule_ids(graph, upstream)
    downstream_ids = _rule_ids(graph, downstream)
    return any(
        edge["u"] in upstream_ids and edge["v"] in downstream_ids
        for edge in graph["edges"]
    )


def _endpoint():
    return "results/annotation_sample/annotation/batch/annotation_results.tsv"


def _singular_selection_endpoint():
    return "results/annotation_sample/annotation/qa/selection.tsv"


def _rule_count(output, rule):
    counts = re.findall(rf"^{re.escape(rule)}\s+([0-9]+)$", output, flags=re.MULTILINE)
    return int(counts[0]) if counts else 0


def _seed_checkpoint(workdir, rows, *, qa=True):
    annotation_root = workdir / "results" / "annotation_sample" / "annotation"
    binning_root = workdir / "results" / "annotation_sample" / "binning" / "domain"
    directories = {
        domain: binning_root / f"{domain}_bins" for domain in ("prok", "euk")
    }
    for directory in directories.values():
        directory.mkdir(parents=True)

    selection_rows = []
    for bin_id, domain, decision in rows:
        fasta = directories[domain] / f"{bin_id}.fa"
        fasta.write_text(">contig\nACGT\n", encoding="utf-8")
        selection_rows.append(
            {
                "sample": bin_id,
                "path": str(fasta.resolve()),
                "domain": domain,
                "completeness": "90" if qa else "",
                "contamination": "1" if qa else "",
                "min_completeness": "50" if qa else "",
                "max_contamination": "10" if qa else "",
                "decision": decision,
                "reason": "fixture",
            }
        )

    if qa:
        checkm2 = annotation_root / "qa" / "prok" / "checkm2" / "quality_report.tsv"
        checkm2.parent.mkdir(parents=True)
        checkm2.write_text(
            "Name\tCompleteness\tContamination\n"
            + "".join(
                f"{bin_id}\t90\t1\n" for bin_id, domain, _ in rows if domain == "prok"
            ),
            encoding="utf-8",
        )
        eukcc_root = annotation_root / "qa" / "euk" / "eukcc"
        eukcc_root.mkdir(parents=True)
        euk_rows = "".join(
            f"{bin_id}\t90\t1\tunknown\n"
            for bin_id, domain, _ in rows
            if domain == "euk"
        )
        (eukcc_root / "eukcc.csv").write_text(
            "bin\tcompleteness\tcontamination\tncbi_lng\n" + euk_rows,
            encoding="utf-8",
        )
        (eukcc_root / "bad_quality.csv").write_text(
            "bin\tcompleteness\tcontamination\tncbi_lng\n", encoding="utf-8"
        )
        (eukcc_root / "missing_marker_genes.txt").write_text("bin\n", encoding="utf-8")

    selection = annotation_root / "qa" / "selection.tsv"
    selection.parent.mkdir(parents=True, exist_ok=True)
    with selection.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SELECTION_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(selection_rows)


def _seed_singular_checkpoint(workdir, domain, *, qa=True):
    annotation_root = workdir / "results" / "annotation_sample" / "annotation"
    assembly = (
        workdir
        / "results"
        / "annotation_sample"
        / "assembly"
        / "assembly-spades"
        / "assembly.fasta"
    )
    selection = annotation_root / "qa" / "selection.tsv"
    selection.parent.mkdir(parents=True, exist_ok=True)
    with selection.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SELECTION_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample": "annotation_sample",
                "path": str(assembly.resolve()),
                "domain": domain,
                "completeness": "90" if qa else "",
                "contamination": "1" if qa else "",
                "min_completeness": "50" if qa else "",
                "max_contamination": "10" if qa else "",
                "decision": "report_only",
                "reason": "fixture",
            }
        )

    if domain == "prok":
        if qa:
            report = annotation_root / "qa" / "prok" / "checkm2" / "quality_report.tsv"
            report.parent.mkdir(parents=True)
            report.write_text(
                "Name\tCompleteness\tContamination\nannotation_sample\t90\t1\n",
                encoding="utf-8",
            )
    else:
        if qa:
            report_root = annotation_root / "qa" / "euk" / "eukcc"
            report_root.mkdir(parents=True)
            (report_root / "eukcc.csv").write_text(
                "bin\tcompleteness\tcontamination\tncbi_lng\n"
                "annotation_sample\t90\t1\tunknown\n",
                encoding="utf-8",
            )
            (report_root / "bad_quality.csv").write_text(
                "bin\tcompleteness\tcontamination\tncbi_lng\n",
                encoding="utf-8",
            )
            (report_root / "missing_marker_genes.txt").write_text(
                "bin\n", encoding="utf-8"
            )


def test_precheckpoint_dag_wires_binning_directly_to_annotation_qa_and_selection(
    tmp_path, all_modules_root
):
    workdir = tmp_path / "precheckpoint"
    workdir.mkdir()
    config_path = _write_config(workdir, all_modules_root)

    graph = _graph(_run(workdir, config_path, "--d3dag", _endpoint()))
    assert _has_edge(
        graph, "binning_separate_bins_by_domain", "annotation_qa_checkm2_prok"
    )
    assert _has_edge(
        graph, "binning_separate_bins_by_domain", "annotation_qa_eukcc_euk"
    )
    assert _has_edge(
        graph, "binning_separate_bins_by_domain", "annotation_qa_selection"
    )
    assert _has_edge(graph, "annotation_qa_checkm2_prok", "annotation_qa_selection")
    assert _has_edge(graph, "annotation_qa_eukcc_euk", "annotation_qa_selection")
    rendered = json.dumps(graph)
    assert "unknown_bins" not in rendered

    listed = _run(workdir, config_path, "--list-rules")
    assert listed.returncode == 0, listed.stdout + listed.stderr
    rules = set(listed.stdout.splitlines())
    assert "annotation_all" not in rules
    assert "metabolic_modeling_complete" not in rules

    harness = workdir / "inspect-context.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("ASSEMBLY_MODULE_CONTEXT=" + assembly.GENOME_CONTEXT)\n'
        'print("ANNOTATION_INPUT_MODE=" + annotation.INPUT_MODE)\n',
        encoding="utf-8",
    )
    inspected = _run(workdir, config_path, "--list-rules", snakefile=harness)
    assert inspected.returncode == 0, inspected.stdout + inspected.stderr
    assert "ASSEMBLY_MODULE_CONTEXT=metagenome" in inspected.stdout
    assert "ANNOTATION_INPUT_MODE=metagenome" in inspected.stdout


@pytest.mark.parametrize(
    ("genome_context", "qa_rule", "annotation_rule"),
    [
        ("prok", "annotation_qa_checkm2_prok", "annotation_protein_summary"),
        ("euk", "annotation_qa_eukcc_euk", "annotation_funannotate2"),
    ],
)
@pytest.mark.parametrize(
    ("qa_enabled", "filter_enabled"),
    [(True, True), (True, False), (False, False)],
    ids=["qa-filter", "qa-report-only", "qa-bypass"],
)
def test_singular_context_routes_assembly_directly_to_annotation_without_binning(
    tmp_path,
    all_modules_root,
    genome_context,
    qa_rule,
    annotation_rule,
    qa_enabled,
    filter_enabled,
):
    workdir = tmp_path / f"{genome_context}-{qa_enabled}-{filter_enabled}"
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        genome_context=genome_context,
        qa=qa_enabled,
        qa_filter=filter_enabled,
        include_binning=False,
    )

    graph = _graph(
        _run(workdir, config_path, "--d3dag", _singular_selection_endpoint())
    )
    rules = {node["value"]["rule"] for node in graph["nodes"]}
    assert not any(rule.startswith("binning_") for rule in rules)
    assert not any("metagenome" in rule for rule in rules)
    if qa_enabled:
        assert _has_edge(graph, "assembly_assemble_spades", qa_rule)
        assert _has_edge(graph, qa_rule, "annotation_qa_selection")
    else:
        assert _has_edge(graph, "assembly_assemble_spades", "annotation_qa_selection")

    harness = workdir / "inspect-main.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("ASSEMBLY_MODULE_CONTEXT=" + assembly.GENOME_CONTEXT)\n'
        'print("ANNOTATION_INPUT_MODE=" + annotation.INPUT_MODE)\n'
        'print("ANNOTATION_DOMAIN=" + annotation.ISOLATE_DOMAIN)\n'
        'print("MAIN_TARGETS=" + "\\n".join(expected_main_outputs()))\n',
        encoding="utf-8",
    )
    listed = _run(workdir, config_path, "--list-rules", snakefile=harness)
    assert listed.returncode == 0, listed.stdout + listed.stderr
    assert "ASSEMBLY_MODULE_CONTEXT=isolate" in listed.stdout
    assert "ANNOTATION_INPUT_MODE=isolate" in listed.stdout
    assert f"ANNOTATION_DOMAIN={genome_context}" in listed.stdout
    assert annotation_rule in listed.stdout.splitlines()
    assert not any(rule.startswith("binning_") for rule in listed.stdout.splitlines())
    targets = listed.stdout.split("MAIN_TARGETS=", 1)[1]
    assert "/assembly/assembly-spades/assembly.fasta" in targets
    assert "/assembly/qc-quast/report.txt" in targets
    assert "/binning/" not in targets
    assert "/assembly/mapping-" not in targets
    if genome_context == "prok":
        assert "/annotation/protein_annotations.tsv" in targets
        assert "/annotation/funannotate2/" not in targets
    else:
        assert "/annotation/funannotate2/summary.json" in targets
        assert "/annotation/protein_annotations.tsv" not in targets

    _seed_singular_checkpoint(workdir, genome_context, qa=qa_enabled)
    dry_run = _run(workdir, config_path, "--dry-run", "--printshellcmds")
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    output = dry_run.stdout + dry_run.stderr
    assert _rule_count(output, annotation_rule) == 1
    assert "binning_" not in output


@pytest.mark.parametrize(
    ("qa_enabled", "filter_enabled"),
    [(True, False), (False, False)],
    ids=["qa-report-only", "qa-bypass"],
)
def test_supported_nondefault_qa_modes_keep_direct_selection_wiring(
    tmp_path,
    all_modules_root,
    qa_enabled,
    filter_enabled,
):
    workdir = tmp_path / f"qa-{qa_enabled}-{filter_enabled}"
    workdir.mkdir()
    config_path = _write_config(
        workdir,
        all_modules_root,
        qa=qa_enabled,
        qa_filter=filter_enabled,
    )

    graph = _graph(_run(workdir, config_path, "--d3dag", _endpoint()))
    assert _has_edge(
        graph, "binning_separate_bins_by_domain", "annotation_qa_selection"
    )
    rules = {node["value"]["rule"] for node in graph["nodes"]}
    for qa_rule in ("annotation_qa_checkm2_prok", "annotation_qa_eukcc_euk"):
        assert (qa_rule in rules) is qa_enabled


@pytest.mark.parametrize(
    ("rows", "expected", "excluded"),
    [
        (
            [("prok_pass", "prok", "pass"), ("euk_pass", "euk", "pass")],
            {
                "annotation_bakta_bin",
                "annotation_recognizer_prok_bin",
                "annotation_upimapi_bin",
                "annotation_protein_summary_bin",
                "annotation_metaeuk_bin",
                "annotation_recognizer_euk_bin",
            },
            set(),
        ),
        (
            [("prok_pass", "prok", "pass")],
            {
                "annotation_bakta_bin",
                "annotation_recognizer_prok_bin",
                "annotation_upimapi_bin",
                "annotation_protein_summary_bin",
            },
            {"annotation_metaeuk_bin", "annotation_recognizer_euk_bin"},
        ),
        (
            [("prok_fail", "prok", "fail"), ("euk_fail", "euk", "fail")],
            set(),
            {
                "annotation_bakta_bin",
                "annotation_recognizer_prok_bin",
                "annotation_upimapi_bin",
                "annotation_protein_summary_bin",
                "annotation_metaeuk_bin",
                "annotation_recognizer_euk_bin",
            },
        ),
        (
            [],
            set(),
            {
                "annotation_bakta_bin",
                "annotation_recognizer_prok_bin",
                "annotation_upimapi_bin",
                "annotation_protein_summary_bin",
                "annotation_metaeuk_bin",
                "annotation_recognizer_euk_bin",
            },
        ),
    ],
    ids=["mixed", "empty-euk-domain", "all-filtered", "empty-domains"],
)
def test_resolved_checkpoint_expands_only_selected_per_bin_routes(
    tmp_path,
    all_modules_root,
    rows,
    expected,
    excluded,
):
    workdir = tmp_path / "resolved"
    workdir.mkdir()
    config_path = _write_config(workdir, all_modules_root)
    _seed_checkpoint(workdir, rows)

    result = _run(workdir, config_path, "--dry-run", "--printshellcmds", _endpoint())
    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    for rule in expected:
        assert _rule_count(output, rule) == 1
    for rule in excluded:
        assert _rule_count(output, rule) == 0
    assert _rule_count(output, "annotation_metagenome_annotation_inventory") == 1
    if "annotation_protein_summary_bin" in expected:
        assert "batch/prok/prok_pass/bakta" in output
        assert "batch/prok/prok_pass/recognizer" in output
        assert "batch/prok/prok_pass/upimapi" in output
        assert "batch/prok/prok_pass/protein_annotations.tsv" in output
    if "annotation_recognizer_euk_bin" in expected:
        assert "batch/euk/euk_pass/metaeuk" in output
        assert "batch/euk/euk_pass/recognizer" in output


@pytest.mark.parametrize("enabled", [False, True], ids=["default-off", "enabled"])
def test_gtdbtk_default_exclusion_and_enabled_post_checkpoint_route(
    tmp_path, all_modules_root, enabled
):
    workdir = tmp_path / f"gtdbtk-{enabled}"
    workdir.mkdir()
    config_path = _write_config(workdir, all_modules_root, gtdbtk=enabled)
    _seed_checkpoint(workdir, [("prok_pass", "prok", "pass")])

    result = _run(workdir, config_path, "--dry-run", "--printshellcmds", _endpoint())
    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    assert (_rule_count(output, "annotation_gtdbtk") == 1) is enabled
    assert "annotation_stage_gtdbtk_genomes" not in output
    assert "annotation/gtdbtk_genomes" not in output
    assert "annotation/gtdbtk_genomes.log" not in output
    if enabled:
        assert "annotation/gtdbtk/classification.tsv" in output
        inventory_job = output.split(
            "rule annotation_metagenome_annotation_inventory:", 1
        )[1].split("Shell command:", 1)[0]
        assert "annotation/gtdbtk/classification.tsv" in inventory_job


@pytest.mark.parametrize(
    ("database", "selected", "database_child", "bakta_type", "bakta_subpath"),
    [
        (None, "UniRef90", "uniref90", None, "bakta_db/db"),
        (
            "  UniProtKB/Swiss-Prot  ",
            "UniProtKB/Swiss-Prot",
            "catalog-UniProtKB%2FSwiss-Prot",
            "light",
            "bakta_db/light/db",
        ),
    ],
    ids=["normal-defaults", "alternate-catalog-and-light"],
)
def test_main_imported_annotation_resolves_database_selections_and_paths(
    tmp_path,
    monkeypatch,
    all_modules_root,
    database,
    selected,
    database_child,
    bakta_type,
    bakta_subpath,
):
    workdir = tmp_path / "database-selection"
    workdir.mkdir()
    config_path = _write_config(workdir, all_modules_root)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["annotation"]["bakta"]["db"] = ""
    config["annotation"]["metaeuk"]["db"] = ""
    config["annotation"]["metaeuk"]["cache_dir"] = ""
    if bakta_type is None:
        config["annotation"]["bakta"].pop("type")
        expected_bakta_type = "full"
    else:
        config["annotation"]["bakta"]["type"] = bakta_type
        expected_bakta_type = bakta_type
    if database is None:
        config["annotation"]["metaeuk"].pop("database")
    else:
        config["annotation"]["metaeuk"]["database"] = database
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    resources = workdir / "selected_resources"
    monkeypatch.delenv("BAKTA_DB", raising=False)
    monkeypatch.delenv("METAEUK_DB", raising=False)
    monkeypatch.setenv("MEEW_RESOURCES", str(resources))
    _seed_checkpoint(
        workdir,
        [("prok_pass", "prok", "pass"), ("euk_pass", "euk", "pass")],
    )
    bakta_database = resources / bakta_subpath
    bakta_database.mkdir(parents=True)
    (bakta_database / "version.json").write_text("{}\n", encoding="utf-8")
    (bakta_database / "bakta.db").write_text("fixture\n", encoding="utf-8")
    metaeuk_prefix = resources / "metaeuk_db" / "catalogs" / database_child / "db"
    metaeuk_prefix.parent.mkdir(parents=True)
    for material in (
        metaeuk_prefix,
        Path(f"{metaeuk_prefix}.index"),
        Path(f"{metaeuk_prefix}.dbtype"),
    ):
        material.write_text("fixture\n", encoding="utf-8")
    (metaeuk_prefix.parent / "catalog-identity.json").write_text(
        "{}\n", encoding="utf-8"
    )

    harness = workdir / "inspect-database-selection.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("METAEUK_DATABASE=" + annotation.METAEUK_DATABASE)\n'
        'print("METAEUK_DB_PREFIX=" + annotation.METAEUK_DB_PREFIX)\n'
        'print("BAKTA_DB_TYPE=" + annotation.BAKTA_DB_TYPE)\n'
        'print("BAKTA_DB_DIR=" + annotation.BAKTA_DB_DIR)\n',
        encoding="utf-8",
    )
    inspected = _run(workdir, config_path, "--list-rules", snakefile=harness)
    assert inspected.returncode == 0, inspected.stdout + inspected.stderr
    assert f"METAEUK_DATABASE={selected}" in inspected.stdout
    assert (
        "METAEUK_DB_PREFIX="
        f"{resources / 'metaeuk_db' / 'catalogs' / database_child / 'db'}"
        in inspected.stdout
    )
    assert f"BAKTA_DB_TYPE={expected_bakta_type}" in inspected.stdout
    assert f"BAKTA_DB_DIR={resources / bakta_subpath}" in inspected.stdout

    dry_run = _run(
        workdir,
        config_path,
        "--dry-run",
        "--printshellcmds",
        _endpoint(),
    )
    assert dry_run.returncode == 0, dry_run.stdout + dry_run.stderr
    rendered = dry_run.stdout + dry_run.stderr
    assert _rule_count(rendered, "annotation_bakta_bin") == 1
    assert _rule_count(rendered, "annotation_metaeuk_bin") == 1
    batch_source = (
        all_modules_root / "MeeW-Annotation" / "workflow" / "rules" / "batch.smk"
    ).read_text(encoding="utf-8")
    for rule, selector in (
        ("bakta_bin", "bakta_database_inputs"),
        ("recognizer_prok_bin", "recognizer_database_inputs"),
        ("upimapi_bin", "upimapi_database_inputs"),
        ("metaeuk_bin", "metaeuk_database_inputs"),
        ("recognizer_euk_bin", "recognizer_database_inputs"),
    ):
        block = batch_source.split(f"rule {rule}:", 1)[1].split("\n\nrule ", 1)[0]
        inputs, outputs = block.split("    output:", 1)
        assert selector in inputs
        assert selector not in outputs


def test_main_forwards_annotation_resource_root_and_preserves_explicit_override(
    tmp_path, monkeypatch, all_modules_root
):
    workdir = tmp_path / "resource-root"
    workdir.mkdir()
    config_path = _write_config(workdir, all_modules_root)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    resources = workdir / "shared-resources"
    explicit_bakta = workdir / "explicit-bakta"
    config["annotation"]["resources_root"] = str(resources)
    config["annotation"]["bakta"]["db"] = str(explicit_bakta)
    config["annotation"]["gtdbtk"]["data_dir"] = ""
    config["annotation"]["funannotate2"]["db_dir"] = ""
    config["annotation"]["metaeuk"]["cache_dir"] = ""
    config["annotation"]["metaeuk"]["db"] = ""
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    for name in (
        "BAKTA_DB",
        "FUNANNOTATE2_DB",
        "GTDBTK_DATA_PATH",
        "METAEUK_DB",
        "MEEW_RESOURCES",
        "XDG_DATA_HOME",
    ):
        monkeypatch.delenv(name, raising=False)

    harness = workdir / "inspect-annotation-resource-root.smk"
    harness.write_text(
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        'print("BAKTA_DB_DIR=" + annotation.BAKTA_DB_DIR)\n'
        'print("GTDBTK_DATA_DIR=" + annotation.GTDBTK_DATA_DIR)\n'
        'print("FUNANNOTATE2_DB_DIR=" + annotation.FUNANNOTATE2_DB_DIR)\n'
        'print("METAEUK_CACHE_DIR=" + annotation.METAEUK_CACHE_DIR)\n'
        'print("METAEUK_DB_PREFIX=" + annotation.METAEUK_DB_PREFIX)\n',
        encoding="utf-8",
    )
    inspected = _run(workdir, config_path, "--list-rules", snakefile=harness)

    assert inspected.returncode == 0, inspected.stdout + inspected.stderr
    assert f"BAKTA_DB_DIR={explicit_bakta}" in inspected.stdout
    assert f"GTDBTK_DATA_DIR={resources / 'gtdbtk_db'}" in inspected.stdout
    assert f"FUNANNOTATE2_DB_DIR={resources / 'funannotate2_db'}" in inspected.stdout
    assert f"METAEUK_CACHE_DIR={resources / 'metaeuk_db'}" in inspected.stdout
    assert (
        "METAEUK_DB_PREFIX="
        f"{resources / 'metaeuk_db' / 'catalogs' / 'uniref90' / 'db'}"
        in inspected.stdout
    )
    assert not resources.exists()
    assert not explicit_bakta.exists()


def test_annotation_schema_accepts_current_modes_and_rejects_main_owned_settings():
    schema = yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema)
    shipped = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    validator.validate(shipped)
    Draft202012Validator(schema["properties"]["annotation"]).validate(
        yaml.safe_load((ROOT / "config" / "config_annotation.yaml").read_text())
    )

    assert shipped["annotation"]["qa"]["enabled"] is True
    assert shipped["annotation"]["qa_filter"] == {
        "enabled": True,
        "min_completeness": 50.0,
        "max_contamination": 10.0,
    }
    assert shipped["annotation"]["gtdbtk"]["enabled"] is False
    assert shipped["annotation"]["bakta"]["type"] == "full"
    assert shipped["annotation"]["metaeuk"]["database"] == "UniRef90"
    assert shipped["annotation"]["resources_root"] == ""
    assert shipped["genome_context"] == "metagenome"
    assert "genome_context" not in shipped["assembly"]

    for qa_enabled, filter_enabled in ((True, True), (True, False), (False, False)):
        configured = copy.deepcopy(shipped)
        configured["annotation"]["qa"]["enabled"] = qa_enabled
        configured["annotation"]["qa_filter"]["enabled"] = filter_enabled
        validator.validate(configured)

    invalid = copy.deepcopy(shipped)
    invalid["annotation"]["qa"]["enabled"] = False
    invalid["annotation"]["qa_filter"]["enabled"] = True
    with pytest.raises(ValidationError):
        validator.validate(invalid)

    omitted_database_selections = copy.deepcopy(shipped)
    omitted_database_selections["annotation"]["bakta"].pop("type")
    omitted_database_selections["annotation"]["metaeuk"].pop("database")
    validator.validate(omitted_database_selections)

    blank_metaeuk_selection = copy.deepcopy(shipped)
    blank_metaeuk_selection["annotation"]["metaeuk"]["database"] = "   "
    validator.validate(blank_metaeuk_selection)

    configured_resources_root = copy.deepcopy(shipped)
    configured_resources_root["annotation"]["resources_root"] = "~/shared/$ROOT"
    validator.validate(configured_resources_root)

    omitted_resources_root = copy.deepcopy(shipped)
    omitted_resources_root["annotation"].pop("resources_root")
    validator.validate(omitted_resources_root)

    invalid_resources_root = copy.deepcopy(shipped)
    invalid_resources_root["annotation"]["resources_root"] = 42
    with pytest.raises(ValidationError):
        validator.validate(invalid_resources_root)

    invalid_bakta_type = copy.deepcopy(shipped)
    invalid_bakta_type["annotation"]["bakta"]["type"] = "compact"
    with pytest.raises(ValidationError):
        validator.validate(invalid_bakta_type)

    for key, value in (
        ("input", {"mode": "metagenome", "sample_sheet": "input.tsv"}),
        ("sample_sheet", "input.tsv"),
        ("stage_sheets", {"enabled": False}),
        ("pathvars", {"results": "other"}),
        ("import_context", {"enabled": True}),
        ("directory_mode", {"enabled": True}),
    ):
        invalid = copy.deepcopy(shipped)
        invalid["annotation"][key] = value
        with pytest.raises(ValidationError):
            validator.validate(invalid)

    for genome_context in ("prok", "euk"):
        singular = copy.deepcopy(shipped)
        singular["genome_context"] = genome_context
        singular.pop("binning")
        validator.validate(singular)

    missing_context = copy.deepcopy(shipped)
    missing_context.pop("genome_context")
    with pytest.raises(ValidationError):
        validator.validate(missing_context)

    invalid_context = copy.deepcopy(shipped)
    invalid_context["genome_context"] = "isolate"
    with pytest.raises(ValidationError):
        validator.validate(invalid_context)

    missing_metagenome_binning = copy.deepcopy(shipped)
    missing_metagenome_binning.pop("binning")
    with pytest.raises(ValidationError):
        validator.validate(missing_metagenome_binning)

    nested_context = copy.deepcopy(shipped)
    nested_context["assembly"]["genome_context"] = "isolate"
    with pytest.raises(ValidationError):
        validator.validate(nested_context)
