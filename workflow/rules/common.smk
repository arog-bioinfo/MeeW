from snakemake.exceptions import WorkflowError

FLYE_LONG_READ_TYPES = preprocessing.FLYE_LONG_READ_TYPES
GENOME_CONTEXT = config["genome_context"]
IS_METAGENOME_CONTEXT = GENOME_CONTEXT == "metagenome"
ASSEMBLY_GENOME_CONTEXT = "metagenome" if IS_METAGENOME_CONTEXT else "isolate"


def validate_main_preprocessing_route():
    rows = tuple(preprocessing.RUN_SHEET_ROWS)
    groups = {row["group"] for row in rows}
    if len(groups) != 1:
        raise WorkflowError(
            "Main MeeW requires exactly one Preprocessing group; found "
            f"{len(groups)} groups: {', '.join(sorted(groups))}."
        )
    if len(rows) > 2:
        raise WorkflowError(
            "Main MeeW accepts one or two normalized Preprocessing rows; "
            f"found {len(rows)} rows for group {rows[0]['group']!r}."
        )

    data_types = [row["data_type"] for row in rows]
    if len(rows) == 1:
        data_type = data_types[0]
        if data_type == "short_pe" or data_type in FLYE_LONG_READ_TYPES:
            return rows
        if data_type == "short_se":
            reason = "Assembly requires paired short reads (data_type=short_pe)"
        elif data_type == "long":
            reason = (
                "legacy data_type=long is not Assembly-compatible because it does "
                "not provide a Flye sequence type"
            )
        else:
            reason = "the data type is not Assembly-compatible"
        raise WorkflowError(
            f"Main MeeW cannot integrate data_type={data_type}: {reason}."
        )

    short_rows = [row for row in rows if row["data_type"] == "short_pe"]
    long_rows = [row for row in rows if row["data_type"] in FLYE_LONG_READ_TYPES]
    if len(short_rows) != 1 or len(long_rows) != 1:
        raise WorkflowError(
            "Main MeeW two-row Preprocessing input requires exactly one "
            "data_type=short_pe row and exactly one Flye-compatible long-read row; "
            f"found data types {data_types!r}."
        )
    return rows


INITIAL_ROWS = validate_main_preprocessing_route()
INITIAL_SAMPLE = INITIAL_ROWS[0]["group"]
INITIAL_IS_HYBRID = len(INITIAL_ROWS) == 2
INITIAL_DATA_TYPE = "hybrid" if INITIAL_IS_HYBRID else INITIAL_ROWS[0]["data_type"]
INITIAL_SHORT_PE_ROW = next(
    (row for row in INITIAL_ROWS if row["data_type"] == "short_pe"), None
)
INITIAL_LONG_ROW = next(
    (row for row in INITIAL_ROWS if row["data_type"] in FLYE_LONG_READ_TYPES), None
)
INITIAL_SHORT_PE_SAMPLE = (
    INITIAL_SHORT_PE_ROW["sample"] if INITIAL_SHORT_PE_ROW else None
)
INITIAL_LONG_SAMPLE = INITIAL_LONG_ROW["sample"] if INITIAL_LONG_ROW else None
INITIAL_LONG_READ_TYPE = INITIAL_LONG_ROW["data_type"] if INITIAL_LONG_ROW else None


def _expand_preprocessing_output(output, **wildcards):
    return [str(path) for path in expand(output, **wildcards)]


def expected_short_processed_reads(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    if INITIAL_IS_HYBRID or INITIAL_DATA_TYPE == "short_pe":
        return _expand_preprocessing_output(
            rules.preprocessing_process_short_pe_fastp.output.trimmed,
            group=sample,
            sample=INITIAL_SHORT_PE_SAMPLE,
        )
    raise WorkflowError(
        f"Sample '{INITIAL_SAMPLE}' does not have paired short reads for Assembly"
    )


def expected_long_processed_reads(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    if INITIAL_IS_HYBRID or INITIAL_DATA_TYPE in FLYE_LONG_READ_TYPES:
        return _expand_preprocessing_output(
            rules.preprocessing_process_long_chopper.output.out_reads,
            group=sample,
            sample=INITIAL_LONG_SAMPLE,
        )[0]
    raise WorkflowError(
        f"Sample '{INITIAL_SAMPLE}' requires a Flye-compatible long-read data_type "
        "for Assembly"
    )


def expected_preprocessing_processed_reads(sample=INITIAL_SAMPLE):
    outputs = []
    if INITIAL_IS_HYBRID or INITIAL_DATA_TYPE == "short_pe":
        outputs.extend(expected_short_processed_reads(sample))
    if INITIAL_IS_HYBRID or INITIAL_DATA_TYPE in FLYE_LONG_READ_TYPES:
        outputs.append(expected_long_processed_reads(sample))
    if outputs:
        return outputs
    raise WorkflowError(f"Unsupported data_type: {INITIAL_DATA_TYPE}")


def expected_preprocessing_qc_outputs(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    outputs = []
    for row in INITIAL_ROWS:
        for qc_stage in ["raw", "processed"]:
            outputs.extend(
                _expand_preprocessing_output(
                    [
                        rules.preprocessing_sequali_qc.output.html,
                        rules.preprocessing_sequali_qc.output.json,
                    ],
                    group=sample,
                    sample=row["sample"],
                    qc_stage=qc_stage,
                )
            )
    return outputs


def expected_preprocessing_outputs(sample=INITIAL_SAMPLE):
    return expected_preprocessing_processed_reads(
        sample
    ) + expected_preprocessing_qc_outputs(sample)


def expected_assembly_type(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    if INITIAL_IS_HYBRID:
        return "hybrid"
    if INITIAL_DATA_TYPE == "short_pe":
        return "short"
    if INITIAL_DATA_TYPE in FLYE_LONG_READ_TYPES:
        return "long"
    raise WorkflowError(
        "Main MeeW Step 2 Assembly wiring requires data_type=short_pe, "
        "a supported two-row hybrid sample, or a Flye-compatible long-read "
        f"data_type ({', '.join(sorted(FLYE_LONG_READ_TYPES))}); "
        f"found data_type={INITIAL_DATA_TYPE}. Legacy data_type=long is "
        "Preprocessing-only and is not Assembly-wirable in main MeeW."
    )


def expected_long_read_type(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    if INITIAL_LONG_READ_TYPE:
        return INITIAL_LONG_READ_TYPE
    raise WorkflowError(
        f"Sample '{INITIAL_SAMPLE}' requires a Flye-compatible long-read data_type "
        "for Assembly"
    )


def expected_hybrid_long_read_key(sample=INITIAL_SAMPLE):
    seq_type = expected_long_read_type(sample)
    if seq_type.startswith("nano-"):
        return "nanopore"
    if seq_type.startswith("pacbio-"):
        return "pacbio"
    raise WorkflowError(
        f"Cannot infer hybrid input key from long_read_type '{seq_type}'"
    )


def _assembly_rule_output(output, sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    return str(output.apply_wildcards({"sample": sample}))


def expected_assembly_output(sample=INITIAL_SAMPLE):
    assembly_type = expected_assembly_type(sample)
    if assembly_type in {"short", "hybrid"}:
        output = workflow.get_rule("assembly_assemble_spades").output.assembly
    elif assembly_type == "long":
        output = workflow.get_rule("assembly_assemble_flye_long").output.assembly
    else:
        raise WorkflowError(f"Unsupported Assembly type: {assembly_type}")
    return _assembly_rule_output(output, sample)


def expected_mapping_bam_output(sample=INITIAL_SAMPLE):
    assembly_type = expected_assembly_type(sample)
    if assembly_type in {"short", "hybrid"}:
        output = workflow.get_rule("assembly_bwamem2_map").output.bam
    elif assembly_type == "long":
        output = workflow.get_rule("assembly_minimap2_map").output.bam
    else:
        raise WorkflowError(f"Unsupported Assembly type: {assembly_type}")
    return _assembly_rule_output(output, sample)


def expected_mapping_bai_output(sample=INITIAL_SAMPLE):
    assembly_type = expected_assembly_type(sample)
    if assembly_type in {"short", "hybrid"}:
        output = workflow.get_rule("assembly_bwamem2_map_index").output.bai
    elif assembly_type == "long":
        output = workflow.get_rule("assembly_minimap2_map_index").output.bai
    else:
        raise WorkflowError(f"Unsupported Assembly type: {assembly_type}")
    return _assembly_rule_output(output, sample)


def expected_mapping_bam_outputs(sample=INITIAL_SAMPLE):
    return [expected_mapping_bam_output(sample)]


def expected_mapping_bai_outputs(sample=INITIAL_SAMPLE):
    return [expected_mapping_bai_output(sample)]


def expected_quast_output(sample=INITIAL_SAMPLE):
    return _assembly_rule_output(
        workflow.get_rule("assembly_quast").output.report_txt, sample
    )


if IS_METAGENOME_CONTEXT:
    BINNING_SELECTED_TOOLS = tuple(config["binning"]["binning"]["tools"])
    BINNING_REFINEMENT_ENABLED = config["binning"]["refinement"]["enabled"]
    BINNING_COVERM_DOMAINS = tuple(
        config["binning"].get("coverm", {}).get("domains", ["prok", "euk"])
    )
else:
    BINNING_SELECTED_TOOLS = ()
    BINNING_REFINEMENT_ENABLED = False
    BINNING_COVERM_DOMAINS = ()


def validate_main_binning_config():
    if not IS_METAGENOME_CONTEXT:
        return
    refinement = config["binning"]["refinement"]
    if (
        refinement["enabled"]
        and refinement["euk_source_tool"] not in BINNING_SELECTED_TOOLS
    ):
        raise WorkflowError(
            "Integrated Binning refinement.euk_source_tool must be one of the "
            f"selected binning.tools; found {refinement['euk_source_tool']!r} with "
            f"{list(BINNING_SELECTED_TOOLS)!r}."
        )


validate_main_binning_config()


def expected_binning_rule_output(
    rule_name, output_name, sample=INITIAL_SAMPLE, **wildcards
):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    output = getattr(workflow.get_rule(rule_name).output, output_name)
    return str(output.apply_wildcards({"sample": sample, **wildcards}))


def expected_refinement_binner_bins_inputs(sample=INITIAL_SAMPLE):
    producer_outputs = {
        "comebin": "binning_comebin",
        "metabinner": "binning_metabinner_bins",
        "semibin2": "binning_semibin2",
        "metabat2": "binning_metabat2",
    }
    return [
        expected_binning_rule_output(producer_outputs[tool], "bins", sample=sample)
        for tool in BINNING_SELECTED_TOOLS
    ]


def expected_refinement_mapping_inputs(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    mappings = workflow.get_rule("binning_prepare_refinement_inputs").output.mappings
    return [str(output.apply_wildcards({"sample": sample})) for output in mappings]


def expected_refinement_assignments(sample=INITIAL_SAMPLE):
    return expected_binning_rule_output(
        "binning_prepare_refinement_inputs",
        "assignments",
        sample=sample,
    )


def expected_domain_bins_dir(domain, sample=INITIAL_SAMPLE):
    if domain not in {"prok", "euk"}:
        raise WorkflowError(
            f"Unsupported integrated CoverM domain '{domain}'; expected prok or euk."
        )
    return expected_binning_rule_output(
        "binning_separate_bins_by_domain",
        f"{domain}_bins",
        sample=sample,
    )


def expected_domain_prok_bins_dir(sample=INITIAL_SAMPLE):
    return expected_binning_rule_output(
        "binning_separate_bins_by_domain", "prok_bins", sample=sample
    )


def expected_domain_euk_bins_dir(sample=INITIAL_SAMPLE):
    return expected_binning_rule_output(
        "binning_separate_bins_by_domain", "euk_bins", sample=sample
    )


def expected_domain_separation_outputs(sample=INITIAL_SAMPLE):
    return [
        expected_binning_rule_output(
            "binning_separate_bins_by_domain", output_name, sample=sample
        )
        for output_name in (
            "prok_bins",
            "euk_bins",
            "unknown_bins",
            "summary",
            "annotation_samples",
            "rename_map",
        )
    ]


def expected_coverm_outputs(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    # CoverM's imported rule currently exposes one genuinely unnamed output.
    output = workflow.get_rule("binning_coverm_genome_coverage").output[0]
    return [
        str(
            output.apply_wildcards(
                {
                    "sample": sample,
                    "domain": domain,
                }
            )
        )
        for domain in BINNING_COVERM_DOMAINS
    ]


def expected_coverm_prok_relative_abundance(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    output = workflow.get_rule("binning_coverm_genome_coverage").output[0]
    return str(output.apply_wildcards({"sample": sample, "domain": "prok"}))


def expected_binning_outputs(sample=INITIAL_SAMPLE):
    if not IS_METAGENOME_CONTEXT:
        return []
    return expected_domain_separation_outputs(sample) + expected_coverm_outputs(sample)


def _expected_annotation_rule_outputs(rule_name, sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    return [
        str(output.apply_wildcards({"sample": sample}))
        for output in workflow.get_rule(rule_name).output
    ]


def expected_annotation_outputs(sample=INITIAL_SAMPLE):
    if GENOME_CONTEXT == "metagenome":
        return _expected_annotation_rule_outputs(
            "annotation_metagenome_annotation_inventory", sample
        )
    if GENOME_CONTEXT == "prok":
        outputs = []
        for rule_name in (
            "annotation_bakta",
            "annotation_recognizer_prok",
            "annotation_upimapi",
            "annotation_protein_summary",
        ):
            outputs.extend(_expected_annotation_rule_outputs(rule_name, sample))
        if config["annotation"].get("gtdbtk", {}).get("enabled", False):
            outputs.extend(
                _expected_annotation_rule_outputs("annotation_gtdbtk", sample)
            )
        return outputs
    if GENOME_CONTEXT == "euk":
        return _expected_annotation_rule_outputs("annotation_funannotate2", sample)
    raise WorkflowError(f"Unsupported genome_context: {GENOME_CONTEXT}")


def annotation_selection_records(wildcards):
    selection = annotation.checkpoints.qa_selection.get(
        sample=INITIAL_SAMPLE
    ).output.selection
    rows = annotation.read_selection(selection)
    seen = set()
    for row in rows:
        key = (row["sample"], row["domain"])
        if key in seen:
            raise WorkflowError(
                "Duplicate integrated Annotation QA selection row for "
                f"{row['domain']} MAG '{row['sample']}' in {selection}."
            )
        seen.add(key)
    return selection, rows


def expected_metabolic_modeling_mags(wildcards):
    if not IS_METAGENOME_CONTEXT:
        return []
    _, rows = annotation_selection_records(wildcards)
    return sorted(
        row["sample"]
        for row in rows
        if row["domain"] == "prok" and row["decision"] == "pass"
    )


def expected_annotation_bakta_directory(bin_id, sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    output = workflow.get_rule("annotation_bakta_bin").output.results
    return str(output.apply_wildcards({"sample": sample, "bin": bin_id}))


def expected_metabolic_modeling_protein_sources(wildcards):
    return {
        bin_id: expected_annotation_bakta_directory(bin_id)
        for bin_id in expected_metabolic_modeling_mags(wildcards)
    }


def expected_annotation_checkm2_report(sample=INITIAL_SAMPLE):
    if sample != INITIAL_SAMPLE:
        raise WorkflowError(f"Unsupported sample: {sample}")
    output = workflow.get_rule("annotation_qa_checkm2_prok").output.report
    return str(output.apply_wildcards({"sample": sample}))


def expected_metabolic_modeling_gtdb_classification(wildcards):
    if not annotation.GTDBTK_ENABLED:
        return []
    output = workflow.get_rule("annotation_gtdbtk").output.classification
    return str(output.apply_wildcards({"sample": INITIAL_SAMPLE}))


def expected_metabolic_modeling_outputs(wildcards):
    if not IS_METAGENOME_CONTEXT or not expected_metabolic_modeling_mags(wildcards):
        return []
    return metabolic_modeling.FINAL_TARGETS(wildcards)


def expected_assembly_outputs():
    return [expected_assembly_output()]


def expected_main_outputs():
    targets = (
        expected_preprocessing_outputs()
        + expected_assembly_outputs()
        + expected_annotation_outputs()
        + [expected_quast_output()]
    )
    if IS_METAGENOME_CONTEXT:
        targets.extend(expected_mapping_bam_outputs())
        targets.extend(expected_mapping_bai_outputs())
        targets.extend(expected_binning_outputs())
    return targets
