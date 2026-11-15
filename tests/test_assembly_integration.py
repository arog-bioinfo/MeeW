import csv
import io
import json
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"

CASES = [
    {
        "mode": "short",
        "sample": "main_short",
        "header": ["group", "sample", "data_type", "read1", "read2"],
        "rows": [
            [
                "main_short",
                "short",
                "short_pe",
                "inputs/short_R1.fastq",
                "inputs/short_R2.fastq",
            ]
        ],
        "assembly_dir": "assembly-spades",
        "mapping_dir": "mapping-bwamem2",
        "rules": {
            "assembly_assemble_spades",
            "assembly_bwamem2_index",
            "assembly_bwamem2_map",
            "assembly_bwamem2_map_index",
            "assembly_quast",
        },
        "opposite_rules": {
            "assembly_assemble_flye_long",
            "assembly_minimap2_map",
            "assembly_minimap2_map_index",
        },
        "preprocessing_rules": {"preprocessing_process_short_pe_fastp"},
        "opposite_preprocessing_rules": {"preprocessing_process_long_chopper"},
    },
    {
        "mode": "long",
        "sample": "main_long",
        "header": ["group", "sample", "data_type", "read1"],
        "rows": [["main_long", "long", "nano-raw", "inputs/long.fastq"]],
        "assembly_dir": "assembly-flye",
        "mapping_dir": "mapping-minimap2",
        "rules": {
            "assembly_assemble_flye_long",
            "assembly_minimap2_map",
            "assembly_minimap2_map_index",
            "assembly_quast",
        },
        "opposite_rules": {
            "assembly_assemble_spades",
            "assembly_bwamem2_index",
            "assembly_bwamem2_map",
            "assembly_bwamem2_map_index",
        },
        "preprocessing_rules": {"preprocessing_process_long_chopper"},
        "opposite_preprocessing_rules": {"preprocessing_process_short_pe_fastp"},
    },
    {
        "mode": "hybrid",
        "sample": "main_hybrid",
        "header": ["group", "sample", "data_type", "read1", "read2"],
        "rows": [
            [
                "main_hybrid",
                "short",
                "short_pe",
                "inputs/short_R1.fastq",
                "inputs/short_R2.fastq",
            ],
            ["main_hybrid", "long", "pacbio-raw", "inputs/long.fastq", ""],
        ],
        "assembly_dir": "assembly-spades",
        "mapping_dir": "mapping-bwamem2",
        "rules": {
            "assembly_assemble_spades",
            "assembly_bwamem2_index",
            "assembly_bwamem2_map",
            "assembly_bwamem2_map_index",
            "assembly_quast",
        },
        "opposite_rules": {
            "assembly_assemble_flye_long",
            "assembly_minimap2_map",
            "assembly_minimap2_map_index",
        },
        "preprocessing_rules": {
            "preprocessing_process_short_pe_fastp",
            "preprocessing_process_long_chopper",
        },
        "opposite_preprocessing_rules": set(),
    },
]


def _write_case(workdir, case, modules_root):
    run_sheet = workdir / "config" / "run.tsv"
    run_sheet.parent.mkdir(parents=True)
    run_sheet.write_text(
        "\t".join(case["header"])
        + "\n"
        + "".join("\t".join(row) + "\n" for row in case["rows"]),
        encoding="utf-8",
    )

    path_columns = [
        index
        for index, column in enumerate(case["header"])
        if column in {"read1", "read2", "contam"}
    ]
    for row in case["rows"]:
        for index in path_columns:
            if row[index]:
                path = workdir / row[index]
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()

    config = yaml.safe_load(
        (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    assert "sample_sheet" not in config["assembly"]
    config["modules_root"] = str(modules_root)
    config["preprocessing"]["run_sheet"] = str(run_sheet)
    config_path = workdir / "config" / "main.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


def _assembly_targets(case):
    root = f"results/{case['sample']}/assembly"
    return [
        f"{root}/{case['assembly_dir']}/assembly.fasta",
        f"{root}/qc-quast/report.txt",
        f"{root}/{case['mapping_dir']}/alignment.bam.bai",
    ]


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["mode"])
def test_main_imported_assembly_route_uses_preprocessing_contract(
    tmp_path, case, all_modules_root
):
    workdir = tmp_path / case["mode"]
    workdir.mkdir()
    config_path = _write_case(workdir, case, all_modules_root)
    targets = _assembly_targets(case)

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
            "--replace-workflow-config",
            "--nolock",
            "--cores",
            "1",
            "--d3dag",
            *targets,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    graph = json.loads(result.stdout)
    rules = {node["value"]["rule"] for node in graph["nodes"]}
    assert case["rules"] <= rules
    assert rules.isdisjoint(case["opposite_rules"])
    assert case["preprocessing_rules"] <= rules
    assert rules.isdisjoint(case["opposite_preprocessing_rules"])
    assert "assembly_all" not in rules
    assert any(rule.startswith("preprocessing_") for rule in rules)

    summary = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
            "--replace-workflow-config",
            "--nolock",
            "--summary",
            *targets,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert summary.returncode == 0, summary.stdout + summary.stderr
    outputs = {
        row["output_file"]
        for row in csv.DictReader(io.StringIO(summary.stdout), delimiter="\t")
    }
    assert set(targets) <= outputs
    opposite_mapping = (
        "mapping-minimap2"
        if case["mapping_dir"] == "mapping-bwamem2"
        else "mapping-bwamem2"
    )
    assert not any(opposite_mapping in output for output in outputs)
