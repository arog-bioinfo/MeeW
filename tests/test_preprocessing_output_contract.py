import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MODULES_ROOT = Path(
    os.environ.get("MEEW_MODULES_ROOT") or REPOSITORY_ROOT.parent
).resolve()
RUN_SHEET_CONFIG_VALUE = "$MAIN_RUN_SHEET"
CONFIG_SCHEMA = REPOSITORY_ROOT / "workflow" / "schemas" / "config.schema.yaml"
CURRENT_INFERENCE_KEYS = [
    "composition_headcrop_enabled",
    "max_headcrop",
    "headcrop_base_imbalance_threshold",
    "tailcrop_mott_quality_cutoff",
    "mott_positive_weight",
    "mott_negative_weight",
    "tailcrop_min_retained_fraction",
    "average_quality",
    "minimum_length",
]
LEGACY_INFERENCE_KEYS = [
    "tailcrop_min_consecutive_good_bases",
    "tailcrop_mean_quality_threshold",
    "tailcrop_min_bad_suffix_span_bases",
    "head_base_imbalance_threshold",
    "quality_crop_mean_threshold",
    "quality_crop_min_bad_suffix_span_bases",
    "quality_crop_min_retained_fraction",
    "quality_crop_low_quality_bucket_boundary",
    "quality_crop_low_quality_fraction_threshold",
]

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
        "data_type": "short_pe",
        "readsets": ["main_short_short"],
        "short_pe_readset": "main_short_short",
        "long_readset": None,
        "long_read_type": None,
        "processed": [
            "reads/main_short_short_R1.fastq",
            "reads/main_short_short_R2.fastq",
        ],
    },
    {
        "mode": "long",
        "sample": "main_long",
        "header": ["group", "sample", "data_type", "read1"],
        "rows": [["main_long", "long", "nano-raw", "inputs/long.fastq"]],
        "data_type": "nano-raw",
        "readsets": ["main_long_long"],
        "short_pe_readset": None,
        "long_readset": "main_long_long",
        "long_read_type": "nano-raw",
        "processed": ["reads/main_long_long.fastq"],
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
            ["main_hybrid", "long", "nano-raw", "inputs/long.fastq", ""],
        ],
        "data_type": "hybrid",
        "readsets": ["main_hybrid_short", "main_hybrid_long"],
        "short_pe_readset": "main_hybrid_short",
        "long_readset": "main_hybrid_long",
        "long_read_type": "nano-raw",
        "processed": [
            "reads/main_hybrid_short_R1.fastq",
            "reads/main_hybrid_short_R2.fastq",
            "reads/main_hybrid_long.fastq",
        ],
    },
]


def _write_run_sheet(workdir, header, rows):
    run_sheet = workdir / "config" / "run.tsv"
    run_sheet.parent.mkdir(parents=True)
    run_sheet.write_text(
        "\t".join(header) + "\n" + "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )
    return run_sheet


def _touch_raw_inputs(workdir, header, rows):
    path_indexes = [
        index
        for index, column in enumerate(header)
        if column in {"read1", "read2", "contam"}
    ]
    for row in rows:
        for index in path_indexes:
            if row[index]:
                path = workdir / row[index]
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()


def _expected_targets(case):
    result_root = f"results/{case['sample']}/preprocessing"
    targets = {f"{result_root}/{path}" for path in case["processed"]}
    targets.update(
        f"{result_root}/qc/{stage}/{readset}.sequali.{extension}"
        for readset in case["readsets"]
        for stage in ["raw", "processed"]
        for extension in ["html", "json"]
    )
    return targets


def _expected_namespace_state(case):
    return {
        "groups": [case["sample"]],
        "rows": [
            {"group": row[0], "sample": row[1], "data_type": row[2]}
            for row in case["rows"]
        ],
    }


def _expected_assembly_state(case):
    result_root = f"results/{case['sample']}/preprocessing"
    state = {
        "assembly_type": {
            "short": "short",
            "long": "long",
            "hybrid": "hybrid",
        }[case["mode"]]
    }
    if case["short_pe_readset"]:
        state["short_reads"] = [
            f"{result_root}/reads/{case['short_pe_readset']}_R1.fastq",
            f"{result_root}/reads/{case['short_pe_readset']}_R2.fastq",
        ]
    if case["long_readset"]:
        state["long_reads"] = f"{result_root}/reads/{case['long_readset']}.fastq"
        state["long_read_type"] = case["long_read_type"]
    if case["mode"] == "hybrid":
        state["hybrid_long_key"] = "nanopore"
    return state


def _write_harness(tmp_path, workdir, *, omit_mott_weights=False):
    with (REPOSITORY_ROOT / "config" / "config.yaml").open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    config["modules_root"] = str(MODULES_ROOT)
    config["preprocessing"]["run_sheet"] = RUN_SHEET_CONFIG_VALUE
    if omit_mott_weights:
        del config["preprocessing"]["auto_preprocess"]["mott_positive_weight"]
        del config["preprocessing"]["auto_preprocess"]["mott_negative_weight"]

    config_path = workdir / "config" / "main.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    harness = tmp_path / "Preprocessing.smk"
    harness.write_text(
        f"""import json
from types import SimpleNamespace

configfile: {json.dumps(str(config_path))}

include: {json.dumps(str(REPOSITORY_ROOT / 'workflow' / 'rules' / 'module_setup.smk'))}

CONFIGURED_RUN_SHEET = config["preprocessing"]["run_sheet"]

module preprocessing:
    pathvars:
        results="results/{{group,[A-Za-z0-9][A-Za-z0-9._-]*}}/preprocessing",
    snakefile:
        str(MODULE_SNAKEFILES["preprocessing"])
    config:
        config["preprocessing"]

use rule * from preprocessing exclude all as preprocessing_*

INFERENCE_RULE_STATE = []
if preprocessing.AUTO_PREPROCESS_ENABLED:
    inference_rule = rules.preprocessing_infer_auto_preprocess
    for row in preprocessing.RUN_SHEET_ROWS:
        if row["data_type"] not in preprocessing.SHORT_READ_TYPES:
            continue
        wildcards = {{"group": row["group"], "sample": row["sample"]}}
        INFERENCE_RULE_STATE.append(
            {{
                "script": str(inference_rule.rule.script),
                "conda": str(inference_rule.rule.conda_env),
                "inputs": {{
                    name: str(getattr(inference_rule.input, name).apply_wildcards(wildcards))
                    for name in ["sequali_json", "config_schema", "sequali_schema"]
                }},
                "output": str(
                    inference_rule.output.auto_params.apply_wildcards(wildcards)
                ),
                "log": str(inference_rule.log[0].apply_wildcards(wildcards)),
                "params": {{
                    "data_type": inference_rule.params.data_type(
                        SimpleNamespace(**wildcards)
                    ),
                    "auto_preprocess_config": dict(
                        inference_rule.params.auto_preprocess_config
                    ),
                }},
            }}
        )
print(
    "PREPROCESSING_INFERENCE_RULE_STATE="
    + json.dumps(INFERENCE_RULE_STATE, sort_keys=True)
)

MODULE_STATE = {{
    "groups": list(preprocessing.GROUPS),
    "rows": [
        {{key: row[key] for key in ["group", "sample", "data_type"]}}
        for row in preprocessing.RUN_SHEET_ROWS
    ],
}}
print("PREPROCESSING_MODULE_STATE=" + json.dumps(MODULE_STATE, sort_keys=True))

include: {json.dumps(str(REPOSITORY_ROOT / 'workflow' / 'rules' / 'common.smk'))}

MAIN_STATE = {{
    "sample": INITIAL_SAMPLE,
    "data_type": INITIAL_DATA_TYPE,
    "is_hybrid": INITIAL_IS_HYBRID,
    "rows": [
        {{key: row[key] for key in ["group", "sample", "data_type"]}}
        for row in INITIAL_ROWS
    ],
    "short_pe_sample": INITIAL_SHORT_PE_SAMPLE,
    "long_sample": INITIAL_LONG_SAMPLE,
    "long_read_type": INITIAL_LONG_READ_TYPE,
    "flye_long_read_types": sorted(FLYE_LONG_READ_TYPES),
    "configured_run_sheet": CONFIGURED_RUN_SHEET,
}}
print("MAIN_NAMESPACE_STATE=" + json.dumps(MAIN_STATE, sort_keys=True))

CONTRACT_TARGETS = expected_preprocessing_outputs()
for contract_target in CONTRACT_TARGETS:
    print("MAIN_PREPROCESSING_TARGET=" + contract_target)

ASSEMBLY_STATE = {{"assembly_type": expected_assembly_type()}}
if INITIAL_IS_HYBRID or INITIAL_DATA_TYPE == "short_pe":
    ASSEMBLY_STATE["short_reads"] = expected_short_processed_reads()
if INITIAL_IS_HYBRID or INITIAL_DATA_TYPE in FLYE_LONG_READ_TYPES:
    ASSEMBLY_STATE["long_reads"] = expected_long_processed_reads()
    ASSEMBLY_STATE["long_read_type"] = expected_long_read_type()
if INITIAL_IS_HYBRID:
    ASSEMBLY_STATE["hybrid_long_key"] = expected_hybrid_long_read_key()
print("MAIN_ASSEMBLY_STATE=" + json.dumps(ASSEMBLY_STATE, sort_keys=True))

rule preprocessing_contract:
    default_target: True
    input:
        CONTRACT_TARGETS,
""",
        encoding="utf-8",
    )
    return harness


def _snakemake_command(harness, workdir):
    snakemake = shutil.which("snakemake")
    assert snakemake, "snakemake must be installed to validate the module contract"
    return [
        snakemake,
        "--snakefile",
        str(harness),
        "--directory",
        str(workdir),
    ]


def _command_environment():
    return {**os.environ, "MAIN_RUN_SHEET": "config/run.tsv"}


def _reported_json(result, prefix):
    payloads = [
        line.removeprefix(prefix)
        for line in result.stdout.splitlines()
        if line.startswith(prefix)
    ]
    assert payloads, result.stdout + result.stderr
    return json.loads(payloads[-1])


def _assert_inference_rule_state(result, case, expected_config):
    states = _reported_json(result, "PREPROCESSING_INFERENCE_RULE_STATE=")
    if case["mode"] == "long":
        assert states == []
        return

    assert len(states) == 1
    state = states[0]
    row = case["rows"][0]
    group, sample, data_type = row[:3]
    result_root = f"results/{group}/preprocessing"
    assert state["script"] == "../scripts/infer_preprocess_params.py"
    assert state["conda"] == "../envs/auto_preprocess.yaml"
    assert set(state["inputs"]) == {
        "sequali_json",
        "config_schema",
        "sequali_schema",
    }
    assert state["inputs"]["sequali_json"] == (
        f"{result_root}/qc/raw/{group}_{sample}.sequali.json"
    )
    assert str(Path(state["inputs"]["config_schema"]).resolve()).endswith(
        "/MeeW-Preprocessing/workflow/schemas/config.schema.yaml"
    )
    assert str(Path(state["inputs"]["sequali_schema"]).resolve()).endswith(
        "/MeeW-Preprocessing/workflow/schemas/sequali_short_read.schema.yaml"
    )
    assert state["output"] == f"{result_root}/reads/{group}_{sample}.auto_params.json"
    assert state["log"] == (
        f"{result_root}/logs/auto_preprocess/{group}_{sample}.auto_params.log"
    )
    assert state["params"] == {
        "data_type": data_type,
        "auto_preprocess_config": expected_config,
    }


@pytest.mark.parametrize(
    "legacy_key",
    ["config", "short_read", "crop_start_position", "min_crop_length"],
)
def test_main_schema_rejects_legacy_auto_preprocess_keys(legacy_key):
    schema = yaml.safe_load(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    Draft202012Validator(schema).validate(config)
    auto_preprocess = config["preprocessing"]["auto_preprocess"]
    if legacy_key == "short_read":
        config["preprocessing"]["auto_preprocess"] = {
            "enabled": auto_preprocess["enabled"],
            "short_read": {
                key: value for key, value in auto_preprocess.items() if key != "enabled"
            },
        }
    else:
        auto_preprocess[legacy_key] = "legacy.yaml"

    with pytest.raises(ValidationError):
        Draft202012Validator(schema).validate(config)


def test_main_active_auto_preprocess_contract_has_nine_scientific_settings():
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )

    assert len(CURRENT_INFERENCE_KEYS) == 9
    assert set(config["preprocessing"]["auto_preprocess"]) == {
        "enabled",
        *CURRENT_INFERENCE_KEYS,
    }
    assert config["preprocessing"]["auto_preprocess"]["mott_positive_weight"] == 1.0
    assert config["preprocessing"]["auto_preprocess"]["mott_negative_weight"] == 1.0


def test_main_schema_allows_mott_weights_to_be_omitted():
    schema = yaml.safe_load(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    del config["preprocessing"]["auto_preprocess"]["mott_positive_weight"]
    del config["preprocessing"]["auto_preprocess"]["mott_negative_weight"]

    Draft202012Validator(schema).validate(config)


@pytest.mark.parametrize("legacy_key", LEGACY_INFERENCE_KEYS)
@pytest.mark.parametrize("disabled", [False, True])
def test_main_schema_rejects_legacy_inference_key_names(legacy_key, disabled):
    schema = yaml.safe_load(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    if disabled:
        config["preprocessing"]["auto_preprocess"] = {"enabled": False}
    config["preprocessing"]["auto_preprocess"][legacy_key] = 20

    with pytest.raises(ValidationError) as error:
        Draft202012Validator(schema).validate(config)
    assert error.value.validator == "additionalProperties"


def test_old_main_config_without_mott_weights_is_accepted_by_imported_module(
    tmp_path,
):
    case = CASES[0]
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    _write_run_sheet(workdir, case["header"], case["rows"])
    _touch_raw_inputs(workdir, case["header"], case["rows"])
    harness = _write_harness(tmp_path, workdir, omit_mott_weights=True)
    written_config = yaml.safe_load(
        (workdir / "config" / "main.yaml").read_text(encoding="utf-8")
    )
    assert (
        "mott_positive_weight" not in written_config["preprocessing"]["auto_preprocess"]
    )
    assert (
        "mott_negative_weight" not in written_config["preprocessing"]["auto_preprocess"]
    )

    result = subprocess.run(
        _snakemake_command(harness, workdir)
        + [
            "--cores",
            "1",
            "--dry-run",
            "--printshellcmds",
            "preprocessing_contract",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_command_environment(),
    )

    assert result.returncode == 0, result.stdout + result.stderr
    expected_config = dict(written_config["preprocessing"]["auto_preprocess"])
    expected_config.update(
        {
            "mott_positive_weight": 1.0,
            "mott_negative_weight": 1.0,
        }
    )
    _assert_inference_rule_state(result, case, expected_config)


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["mode"])
def test_namespace_drives_main_preprocessing_and_assembly_contract(tmp_path, case):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    _write_run_sheet(workdir, case["header"], case["rows"])
    _touch_raw_inputs(workdir, case["header"], case["rows"])
    harness = _write_harness(tmp_path, workdir)
    base_command = _snakemake_command(harness, workdir)

    result = subprocess.run(
        base_command
        + [
            "--cores",
            "1",
            "--dry-run",
            "--printshellcmds",
            "preprocessing_contract",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_command_environment(),
    )
    assert result.returncode == 0, result.stdout + result.stderr
    main_config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    _assert_inference_rule_state(
        result,
        case,
        main_config["preprocessing"]["auto_preprocess"],
    )

    module_state = _reported_json(result, "PREPROCESSING_MODULE_STATE=")
    assert module_state == _expected_namespace_state(case)
    main_state = _reported_json(result, "MAIN_NAMESPACE_STATE=")
    assert main_state.pop("configured_run_sheet") == RUN_SHEET_CONFIG_VALUE
    assert main_state == {
        "sample": case["sample"],
        "data_type": case["data_type"],
        "is_hybrid": case["mode"] == "hybrid",
        "rows": module_state["rows"],
        "short_pe_sample": (
            case["rows"][0][1] if case["mode"] in {"short", "hybrid"} else None
        ),
        "long_sample": (
            case["rows"][-1][1] if case["mode"] in {"long", "hybrid"} else None
        ),
        "long_read_type": case["long_read_type"],
        "flye_long_read_types": main_state["flye_long_read_types"],
    }
    assert _reported_json(result, "MAIN_ASSEMBLY_STATE=") == (
        _expected_assembly_state(case)
    )

    reported_targets = {
        line.removeprefix("MAIN_PREPROCESSING_TARGET=")
        for line in result.stdout.splitlines()
        if line.startswith("MAIN_PREPROCESSING_TARGET=")
    }
    assert reported_targets == _expected_targets(case)
    assert all("processed_reads" not in target for target in reported_targets)

    listed_rules = subprocess.run(
        base_command + ["--list-rules"],
        check=False,
        capture_output=True,
        text=True,
        env=_command_environment(),
    )
    assert listed_rules.returncode == 0, listed_rules.stdout + listed_rules.stderr
    assert "preprocessing_all" not in listed_rules.stdout.splitlines()


@pytest.mark.parametrize("data_type", ["short_se", "long"])
def test_main_rejects_standalone_only_route_after_module_validation(
    tmp_path, data_type
):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    header = ["group", "sample", "data_type", "read1"]
    rows = [["standalone_only", "reads", data_type, "inputs/reads.fastq"]]
    _write_run_sheet(workdir, header, rows)
    _touch_raw_inputs(workdir, header, rows)
    harness = _write_harness(tmp_path, workdir)

    result = subprocess.run(
        _snakemake_command(harness, workdir) + ["--cores", "1", "--dry-run"],
        check=False,
        capture_output=True,
        text=True,
        env=_command_environment(),
    )

    assert result.returncode != 0
    combined_output = result.stdout + result.stderr
    assert f"data_type={data_type}" in combined_output
    assert "Assembly" in combined_output
    assert "Main MeeW cannot integrate" in combined_output


def test_main_uses_generated_row_sample_to_select_preprocessing_outputs(tmp_path):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    header = ["group", "sample", "data_type", "read1", "read2"]
    rows = [
        [
            "generated_group",
            "",
            "short_pe",
            "inputs/short_R1.fastq",
            "inputs/short_R2.fastq",
        ]
    ]
    _write_run_sheet(workdir, header, rows)
    _touch_raw_inputs(workdir, header, rows)
    harness = _write_harness(tmp_path, workdir)

    result = subprocess.run(
        _snakemake_command(harness, workdir)
        + ["--cores", "1", "--dry-run", "preprocessing_contract"],
        check=False,
        capture_output=True,
        text=True,
        env=_command_environment(),
    )

    assert result.returncode == 0, result.stdout + result.stderr
    module_state = _reported_json(result, "PREPROCESSING_MODULE_STATE=")
    assert module_state["rows"] == [
        {
            "group": "generated_group",
            "sample": "short_pe_1",
            "data_type": "short_pe",
        }
    ]
    assembly_state = _reported_json(result, "MAIN_ASSEMBLY_STATE=")
    assert assembly_state["short_reads"] == [
        "results/generated_group/preprocessing/reads/"
        "generated_group_short_pe_1_R1.fastq",
        "results/generated_group/preprocessing/reads/"
        "generated_group_short_pe_1_R2.fastq",
    ]


@pytest.mark.parametrize(
    ("case_name", "header", "rows", "diagnostic"),
    [
        (
            "multiple-groups",
            ["group", "sample", "data_type", "read1"],
            [
                ["group1", "long", "nano-raw", "inputs/one.fastq"],
                ["group2", "long", "pacbio-raw", "inputs/two.fastq"],
            ],
            "requires exactly one Preprocessing group",
        ),
        (
            "more-than-two-rows",
            ["group", "sample", "data_type", "read1", "read2"],
            [
                [
                    "group1",
                    "short",
                    "short_pe",
                    "inputs/short_R1.fastq",
                    "inputs/short_R2.fastq",
                ],
                ["group1", "long1", "nano-raw", "inputs/long1.fastq", ""],
                ["group1", "long2", "pacbio-raw", "inputs/long2.fastq", ""],
            ],
            "accepts one or two normalized Preprocessing rows",
        ),
        (
            "duplicate-short",
            ["group", "sample", "data_type", "read1", "read2"],
            [
                [
                    "group1",
                    "short1",
                    "short_pe",
                    "inputs/one_R1.fastq",
                    "inputs/one_R2.fastq",
                ],
                [
                    "group1",
                    "short2",
                    "short_pe",
                    "inputs/two_R1.fastq",
                    "inputs/two_R2.fastq",
                ],
            ],
            "exactly one data_type=short_pe row",
        ),
        (
            "duplicate-long",
            ["group", "sample", "data_type", "read1"],
            [
                ["group1", "long1", "nano-raw", "inputs/one.fastq"],
                ["group1", "long2", "pacbio-hifi", "inputs/two.fastq"],
            ],
            "exactly one Flye-compatible long-read row",
        ),
        (
            "short-se-plus-long",
            ["group", "sample", "data_type", "read1"],
            [
                ["group1", "short", "short_se", "inputs/short.fastq"],
                ["group1", "long", "nano-hq", "inputs/long.fastq"],
            ],
            "exactly one data_type=short_pe row",
        ),
        (
            "short-plus-legacy-long",
            ["group", "sample", "data_type", "read1", "read2"],
            [
                [
                    "group1",
                    "short",
                    "short_pe",
                    "inputs/short_R1.fastq",
                    "inputs/short_R2.fastq",
                ],
                ["group1", "long", "long", "inputs/long.fastq", ""],
            ],
            "exactly one Flye-compatible long-read row",
        ),
        (
            "two-standalone-only-modes",
            ["group", "sample", "data_type", "read1"],
            [
                ["group1", "single", "short_se", "inputs/short.fastq"],
                ["group1", "legacy", "long", "inputs/long.fastq"],
            ],
            "two-row Preprocessing input requires exactly one",
        ),
    ],
)
def test_main_rejects_invalid_group_cardinality_and_mode_pairs(
    tmp_path, case_name, header, rows, diagnostic
):
    workdir = tmp_path / case_name
    workdir.mkdir()
    _write_run_sheet(workdir, header, rows)
    _touch_raw_inputs(workdir, header, rows)
    harness = _write_harness(tmp_path, workdir)

    result = subprocess.run(
        _snakemake_command(harness, workdir) + ["--cores", "1", "--dry-run"],
        check=False,
        capture_output=True,
        text=True,
        env=_command_environment(),
    )

    assert result.returncode != 0
    assert diagnostic in result.stdout + result.stderr
