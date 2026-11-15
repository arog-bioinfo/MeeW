## Workflow overview

The main MeeW workflow is an import-based Snakemake DAG. It imports sibling stage sources with group-derived, stage-separated
result directories such as `results/{sample}/assembly` and connects producer
outputs directly to consumers. The required top-level
`genome_context` is the only routing source: `metagenome` (the default) runs
Preprocessing → Assembly → Binning → Annotation and continues through Metabolic
Modeling when Annotation selects passing prokaryotic MAGs, while `prok` and
`euk` run Preprocessing → Assembly → Annotation only.

## Portable paths and module discovery

All module input/config paths expand `~` and environment variables. Relative
user input paths intentionally remain relative to the Snakemake working
directory. Main MeeW discovers the five local module repositories in this
order:

1. `modules_root` in `config/config.yaml`;
2. `MEEW_MODULES_ROOT`;
3. the parent directory containing this repository and side-by-side `MeeW-*`
   repositories.

The selected root must contain `MeeW-Preprocessing`, `MeeW-Assembly`,
`MeeW-Binning`, `MeeW-Annotation`, and `MeeW-Metabolic_Modeling`. This keeps
local sibling development and does not fetch remote workflow modules.
An empty `preprocessing.porechop_abi.ab_initio_config` value uses the
corresponding file from the resolved `MeeW-Preprocessing` repository; an
explicit path still takes precedence. Short-read auto-preprocess settings live
directly under `preprocessing.auto_preprocess` in the single active Main config
as `preprocessing.auto_preprocess.<setting>`. Main passes that section to
MeeW-Preprocessing, which owns its detailed scientific validation.

For each observed mean-quality anchor, Preprocessing computes
`difference = tailcrop_mott_quality_cutoff - mean_quality`. The Mott score
contribution is `difference * mott_positive_weight` when the difference is
positive, `difference * mott_negative_weight` when it is negative, and exactly
`0` when it is zero. Both weights are optional nonnegative numbers; omission
defaults to `1.0`/`1.0`, and the canonical Main config shows those values
explicitly. A zero positive weight removes poor-quality evidence, a zero
negative weight removes good-quality offset, and both zero produce no positive
candidate. The relative ratio is operationally meaningful, while multiplying
both weights by the same positive factor only rescales the cumulative score and
does not change the selected candidate.

Integrated cross-repository CI is intentionally deferred. Coordinated sibling
checkout authentication and immutable revision pinning need a separate release
policy before CI can validate this local sibling layout; the current workflow
must not be treated as that policy.

Top-level `resources_root` is an optional Main-wide persistent resource root.
When nonblank, Main forwards it to Binning, Annotation, and Metabolic Modeling
only when that module's own `resources_root` is blank or omitted. A nonblank
module root overrides Main. Explicit tool paths continue to override either
root. Blank Main and module roots preserve every module's existing ambient
defaults instead of redirecting them through config merging.

Binning selects CheckM 1 in this order: exact nonblank
`binning.checkm1.database_path`, nonblank `binning.resources_root`, forwarded
nonblank Main `resources_root`, `MEEW_RESOURCES`, `$XDG_DATA_HOME/meew`, then
`~/resources`; every root-derived path adds `checkm1_db`. Automatic setup is
enabled by default. Before COMEBin augmentation/training or MetaBinner component
binning/scoring and output cleanup, Binning reuses a complete selected cache or
safely invokes that tool environment's installed CheckM Bioconda post-link hook
in staging, validates bounded native scoring and tool-specific marker readiness,
publishes it atomically, and verifies native registration. A canonical lock in
the concrete environment is retained through the complete CheckM-consuming job;
this supplies one shared lock inode to path aliases visible through the same
filesystem but does not claim local tests prove distributed lock behavior.

Annotation retains its existing explicit-tool/module-root/tool-variable/ambient
precedence and derived resource subdirectories. Metabolic Modeling now accepts
the same optional module root for its automatic CheckM2 path; explicit CheckM2
config wins, a selected module/Main root precedes `CHECKM2DB`, and a blank root
preserves `CHECKM2DB`, `MEEW_RESOURCES`, `$XDG_DATA_HOME/meew`, then
`~/.local/share/meew`. Main still disables and excludes the module's internal
CheckM2 execution. Configured paths expand environment variables and `~`;
relative paths use Snakemake's execution working directory.

CheckM's Bioconda package normally invokes its post-link downloader only while
the environment is installed. Setting `CHECKM_DATA_DIR` on an already cached
environment does not rerun that install hook; `CHECKM_DATA_PATH` is interpreted
by the hook as a local archive, not a database-directory override. Binning's
job-time setup is the distinct runtime mechanism for COMEBin and MetaBinner and
removes ambient `CHECKM_DATA_PATH` so it cannot redirect the configured
destination. Other module database setup behavior remains owned by those
modules, and MeeW does not bundle databases.

The older `config/config_preprocessing.yaml`, `config/config_assembly.yaml`,
`config/config_binning.yaml`, `config/config_annotation.yaml`, and
`config/config_metabolic_modeling.yaml` fragments are retained only as
deprecated stage-specific examples. They are not loaded by the active
orchestrator; use `config/config.yaml`.

Step 2 wires Preprocessing into Assembly for the initial sample. Main derives
the Assembly module setting from the top-level context: `metagenome` maps to
Assembly `genome_context: metagenome`, while `prok` and `euk` both map to
Assembly `genome_context: isolate`. Do not add `assembly.genome_context`; Main's
schema rejects that conflicting second source.
Main passes the configured `preprocessing.run_sheet` string into MeeW-Preprocessing unchanged. MeeW-Preprocessing alone expands that value, parses the TSV, validates and normalizes rows, generates blank row sample IDs, and exports normalized state. Each row requires the headers `group`, `sample`, `data_type`, and `read1`; `group` is nonblank, while sample cells may be blank. `read2` and `contam` are valid only where the Preprocessing schema permits them. Legacy identifier columns (`dataset`, `name`, `id`, and `readset`) are not accepted by that schema.
After module evaluation and rule import, Main consumes `preprocessing.RUN_SHEET_ROWS`, `preprocessing.FLYE_LONG_READ_TYPES`, and imported output objects without reparsing the run sheet. Main requires exactly one group with one or two rows. One row must be `short_pe` or a Flye-compatible long type (`nano-raw`, `nano-corr`, `nano-hq`, `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`). Two rows must be exactly one `short_pe` plus one Flye-compatible long row. Multiple groups, more than two rows, `short_se`, legacy `long`, duplicate modes, and all other pairs are rejected. Row sample values select Preprocessing outputs; the group becomes the downstream sample and result identity.
Main derives processed-read and Sequali targets from imported Preprocessing rule outputs. Every row uses the exact `{group}_{row_sample}` basename beneath `results/{group}/preprocessing/`; outputs produced alongside final FASTQs and Sequali HTML/JSON targets, such as logs and fastp reports, are not separately forced.
Main uses Snakemake's module validation control for Assembly's standalone-only required configuration and selectively imports only the current route's rules. It directly overrides sample-dependent inputs and parameters from `preprocessing.*`, so Assembly's lazily loaded standalone record is not requested and Main does not create or parse an Assembly input sheet. Short and hybrid use SPAdes plus BWA-MEM2, while long uses Flye plus minimap2. The minimap2 preset comes from the normalized long row's `data_type`; `assembly.minimap2_extra` must not set `-x`.
For `genome_context: metagenome`, Step 3 wires the selected current Assembly outputs directly into mandatory Binning. Short/hybrid contigs are `results/{sample}/assembly/assembly-spades/assembly.fasta` with the exact `mapping-bwamem2/alignment.bam`; long contigs are `results/{sample}/assembly/assembly-flye/assembly.fasta` with the exact `mapping-minimap2/alignment.bam`. Main passes no Binning `sample_sheet`. Singular `prok`/`euk` contexts do not import Binning or force mapping outputs; they route the selected Assembly FASTA directly to Annotation. QUAST remains `results/{sample}/assembly/qc-quast/report.txt` in every context.

With `binning.refinement.enabled: false`, `binning.binning.tools` must contain exactly one of `comebin`, `metabinner`, `semibin2`, or `metabat2`. With refinement enabled, use `method: das_tool`, select at least two distinct tools, and set `euk_source_tool` to one selected tool. Each selected binner is reached through MeeW-Binning's imported internal producers rather than a Main reimplementation. Main injects SemiBin2 `sequencing_type`: a short Assembly route uses `short_read`; long and hybrid routes use `long_read`. Users cannot set a conflicting value.

Main preserves the imported COMEBin and MetaBinner rules' module-owned
runtime-helper source inputs while overriding only Assembly contigs/BAM and
MetaBinner profile material. Dry runs therefore show the selected CheckM 1
destination and helper sources but perform no database setup, registration,
download, or filesystem initialization.

Final Binning targets come from the imported domain-separation output object: `prok_bins`, `euk_bins`, `unknown_bins`, `domain_assignment.tsv`, `annotation_samples.tsv`, and `bin_rename_map.tsv`. Unknown bins are retained but excluded from downstream coverage and Annotation. Integrated CoverM is mandatory, supports only `input_mode: bam`, consumes the exact Assembly BAM, and always targets the configured unique nonempty subset of `prok`/`euk`; the removed `binning.coverm.enabled` switch, raw-read input fields, user-supplied BAM paths, and unknown-domain coverage are rejected.

MAG QA is not a Binning responsibility. For metagenomes, Main injects the initial sample and exact future imported Binning `prok_bins`/`euk_bins` outputs into Annotation's private module context. For singular runs, it injects the sample, `prok`/`euk` domain, and future selected Assembly FASTA. Future producer outputs are declared without parse-time existence checks. Main imports only the context-appropriate Annotation rules. Users cannot set Annotation `input`, `sample_sheet`, `import_context`, legacy `directory_mode`, or `pathvars` values in Main config. Annotation's strict standalone one-record samplesheet contracts are unchanged.

The shipped `annotation` defaults are `qa.enabled: true`, `qa_filter.enabled: true`, `min_completeness: 50.0`, `max_contamination: 10.0`, and `gtdbtk.enabled: false`. QA on/filter off retains all bins after reporting, and QA off/filter off bypasses scientific QA while retaining the selection contract. QA off/filter on is invalid. Singular QA is always report-only. Passing metagenome prokaryotic bins receive Bakta, reCOGnizer, UPIMAPI, and protein summaries; passing eukaryotic bins receive MetaEuk and reCOGnizer. Singular `prok` uses the prokaryotic toolchain, while singular `euk` uses Funannotate2 and requires a valid species-training, params-JSON, or pretrained-model configuration. Enabled GTDB-Tk receives passing prokaryotic material only.

For integrated Metabolic Modeling, enabled GTDB-Tk's named consolidated `classification.tsv` output is also a declared input to every selected CarveMe job. Main only wires that imported material; the Metabolic Modeling module resolves exact `user_genome` rows at job time and passes domain `bacteria` or `archaea` as the same-spelled CarveMe `--universe`. Blank/unknown domains or no exact row visibly fall back to `bacteria`; corrupt, duplicate, missing, or unreadable promised classification input fails. GTDB-Tk is not auto-enabled, and its disabled route adds no classification dependency and uses the bacterial default.

Annotation database selection is distinct from location. `annotation.bakta.type`
is optionally `full` (the normal default) or `light`; `annotation.bakta.db` is
an optional explicit local directory. Without a higher-precedence location, the
respective resource-root subpaths are `bakta_db/db` and `bakta_db/light/db`.
`annotation.metaeuk.database` is an optional free-form
MMseqs2 catalog name: omission, empty text, or whitespace uses `UniRef90`, and
a custom value is trimmed and otherwise passed unchanged as one argument.
`annotation.metaeuk.db` remains an optional explicit local prefix, while
`annotation.metaeuk.cache_dir` is the cache root. A blank prefix derives
`catalogs/uniref90/db` for `UniRef90`; another source derives
`catalogs/catalog-<urllib.parse.quote(trimmed_source, safe="")>/db` beneath that
cache, such as `catalogs/catalog-UniProtKB%2FSwiss-Prot/db`. Annotation database
payloads are persistent adapter-managed caches, not Snakemake outputs or setup
targets. Ready matching caches are reused without installer/network calls on
annotation reruns. Bakta validates schema/type/native material. Managed MetaEuk
leaves validate the complete data/index/dbtype family and a catalog identity
record. Complete explicit or historical flat prefixes remain in place; absent
identity is reported as a user assertion, known mismatches fail, and missing
external prefixes are not populated automatically.

`annotation.funannotate2.databases` retains its existing nonempty, unique
YAML-list behavior. The shipped `- all` entry (equivalent to `[all]`) installs
all databases. A named subset such as the following is accepted and expanded
to separate argv elements:

```yaml
databases:
  - uniprot
  - pfam
  - dbCAN
```

That list becomes `funannotate2 install -d uniprot pfam dbCAN`. Names are
case-sensitive, especially `dbCAN`; allowed values are `all`, `merops`,
`uniprot`, `dbCAN`, `pfam`, `go`, `mibig`, `interpro`, `gene2product`, and
`mito`. Upstream accepts `all` mixed with other names, and MeeW intentionally
does not reject those lists or change upstream behavior. In Funannotate2
26.6.21, the installer checks membership of `all` and replaces the entire
requested list with all named databases, making extra entries redundant rather
than additive.

The Annotation endpoint is the imported output `results/{sample}/annotation/batch/annotation_results.tsv`. Its DAG dependencies include Annotation's QA selection and every selected per-bin material. Empty domain directories, domains with no bins, and all-filtered selections remain valid and produce no per-bin route for that set. Singular final targets are the imported prokaryotic annotation materials or Funannotate2 outputs, with no Binning outputs. Focused integration tests use all five real sibling module sources.

For metagenomes, Main reads the resolved imported Annotation QA checkpoint with `annotation.read_selection` and selects only rows whose domain is `prok` and decision is `pass`. It formats each selected `annotation_bakta_bin` named `results` output with the real `sample` and `bin` wildcards after pathvar relocation. Those exact directories are material inputs of both the imported Metabolic Modeling external-QA checkpoint and corresponding CarveMe jobs. At job time, module-owned code validates the fixed `annotation.faa` child; Main does not glob, discover, stage, or infer MAG IDs from files.

Main clears standalone Metabolic Modeling FAA/sample-sheet/glob/directory/Bakta-discovery inputs, injects the initial sample, disables and excludes its internal CheckM2 rules, and relocates outputs to `results/{sample}/metabolic_modeling`. When Annotation QA and `metabolic_modeling.qa.enabled` are both true, Main directly supplies Annotation's `qa_checkm2_prok` report with canonical `Name`, `Completeness`, and `Contamination` columns. Annotation QA off forces integrated Metabolic Modeling QA off. The Main schema accepts solver (CPLEX remains the default), reconstruction media, execution mock, QA enable/thresholds, Memote, community tool selection, the shared `simulation_medium` selection, SMETANA, and SteadyCom settings, while rejecting Main-owned adapter paths, identity, QA columns, pathvars, internal CheckM2, and user-supplied `metabolic_modeling.micom.abundance_tsv`. Main forwards `simulation_medium` and relies on module-owned validation and implementation rather than parsing its scientific contents.

Main always injects the exact imported prokaryotic CoverM output `results/{sample}/binning/coverm/prok_genome_coverage.tsv` as MICOM's declared abundance input. The module owns CoverM-table parsing and retained-model normalization; Main performs no scientific parsing. Enabled metagenome community MICOM therefore requires `prok` in `binning.coverm.domains` and `relative_abundance` in `binning.coverm.methods`. Singular routes and configurations with community or MICOM disabled do not require that pair.

The default `metabolic_modeling.simulation_medium` uses `source: discover`, absolute solo `min_growth: 0.1`, existing-exchange opening at `exchange_bound: 1000`, and unweighted `aggregation: sum`. `aggregation: max` is the other discovery choice. `source: predefined` schedules no individual medium optimization and selects either bundled `M9` (20 nutrients, including oxygen, each uptake 100) or a user `nutrient_id,uptake` CSV path. M9 is an aerobic computational flux-budget preset, not an experimentally measured recipe. This selection affects MICOM and SteadyCom only; SMETANA remains unchanged.

For at least one qualifying prokaryotic MAG, checkpoint-aware `rule all` uses the imported module namespace's existing concrete `FINAL_TARGETS` through external QA, CarveMe, Memote, filtering, and configured community rules. Euk-only, empty, and all-filtered metagenomes add no Metabolic Modeling targets and finish at Annotation. Singular `prok` and `euk` routes do not evaluate or import Metabolic Modeling.

## Standalone module execution

Run the module repositories one at a time in order:

1. `MeeW-Preprocessing`
2. `MeeW-Assembly`
3. `MeeW-Binning`
4. `MeeW-Annotation`
5. `MeeW-Metabolic_Modeling`

Each standalone repository retains its user-authored input-sheet mechanism. To chain standalone runs, author the next module's input sheet from the prior module's documented scientific outputs. Main instead targets stage-owned outputs directly under `results/{sample}/{stage}`. Binning and conditional direct Metabolic Modeling are metagenome-only.

At a high level:

- Preprocessing produces row-scoped FASTQs and QC reports; configure Assembly with the selected read outputs.
- Assembly produces the route-selected assembly and mapping outputs; configure Binning from those outputs.
- Binning separates prokaryotic/eukaryotic/unknown bins; configure Annotation from the final domain directories. QA is Annotation-owned.
- Metabolic Modeling retains its generic user-authored `input.faa_samples` TSV and other direct FAA discovery modes.

## Main MeeW role

Main MeeW coordinates context-selected module rules through direct imported output contracts. Its default metagenome route conditionally continues through directly imported Metabolic Modeling rules while preserving standalone user-input mechanisms.
