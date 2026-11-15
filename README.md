# Snakemake pipeline: MeeW

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A59.12.0-039475?logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/conda/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<!-- TEMPORARY CATALOG OPT-OUT: later replace every standalone "pipeline" with the standard term formed by "work" + "flow", then remove this comment. -->

Main MeeW is a Snakemake orchestrator for metagenome and singular prokaryotic or eukaryotic genome analysis. It connects five independently runnable local pipeline modules while leaving their scientific implementations module-owned.

> [!WARNING]
> **Research preview:** MeeW is under active development and is not a deployment-ready release. Validate its results and resource requirements for your use case.

## Usage

The independently runnable modules are [Preprocessing](https://github.com/arog-bioinfo/MeeW-Preprocessing), [Assembly](https://github.com/arog-bioinfo/MeeW-Assembly), [Binning](https://github.com/arog-bioinfo/MeeW-Binning), [Annotation](https://github.com/arog-bioinfo/MeeW-Annotation), and [Metabolic Modeling](https://github.com/arog-bioinfo/MeeW-Metabolic_Modeling).

Each Main run accepts one group containing one Assembly-compatible Preprocessing row (`short_pe` or a Flye-compatible long-read type), or two rows containing exactly one of each for a hybrid assembly. The `metagenome` route runs Preprocessing → Assembly → Binning → Annotation and continues to Metabolic Modeling only for passing prokaryotic MAGs selected by Annotation. Metagenomes with no qualifying MAGs finish at Annotation. Singular `prok` and `euk` routes skip Binning and Metabolic Modeling and finish at Annotation. MultiQC remains planned and is not an implemented stage.

Install Snakemake `>=9.12.0` and Conda. From the repository root, bootstrap missing local modules or verify existing checkouts with:

```bash
bash scripts/setup_modules.sh
```

The helper clones a missing module's mutable `main` branch and is not an immutable release installation; existing valid module checkouts are only verified and retained. Local module sources also do not make the full pipeline offline.

Edit [`config/config.yaml`](config/config.yaml) for the route, tools, and resources, then replace the example paths in [`config/preprocessing_run_sheet.tsv`](config/preprocessing_run_sheet.tsv). Databases are not bundled, and resource requirements depend on the selected tools. Native MICOM and SteadyCom execution requires a separately eligible, installed, and licensed CPLEX installation.

Results are organized by group and stage below `results/`. See the [configuration guide](config/README.md) and the linked module documentation for detailed inputs, algorithms, resources, and output contracts.

## Deployment options

Conda is the intended software deployment method. After configuration, inspect the planned jobs and then execute the pipeline from the repository root with:

```bash
snakemake --cores 2 --dry-run
snakemake \
  --cores 2 \
  --software-deployment-method conda
```

Adjust cores and resource overrides for the selected tools, databases, and execution environment.

## Execution profiles

Main supplies no general configured execution profile. Follow the [profile guidance](profiles/README.md) and add or select a user-maintained profile for executor and site-specific settings.

## License

The pipeline is distributed under the [MIT License](LICENSE).

## References

No MeeW-specific DOI is declared. Cite the repository version used and the selected module repositories, scientific tools, solvers, and databases, together with Snakemake:

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
