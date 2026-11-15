import sys
from pathlib import Path

from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

_MAIN_PATHS_SOURCE = Path(workflow.source_path("../scripts/meew_main_paths.py"))
sys.path.insert(0, str(_MAIN_PATHS_SOURCE.parent))
from meew_main_paths import expand_path as expand_configured_path
from meew_main_paths import resolve_modules_root


def _resolve_modules_root():
    try:
        return resolve_modules_root(
            config.get("modules_root", ""), workflow.basedir, cwd=Path.cwd()
        )
    except FileNotFoundError as error:
        raise WorkflowError(
            f"{error}. Set config.modules_root or "
            "MEEW_MODULES_ROOT to the directory containing the side-by-side module "
            "repositories."
        ) from error


def validate_integrated_micom_coverm_config():
    if config.get("genome_context") != "metagenome":
        return
    metabolic_modeling = config.get("metabolic_modeling", {})
    if not isinstance(metabolic_modeling, dict):
        return
    community = metabolic_modeling.get("community", {})
    if not isinstance(community, dict):
        return
    tools = community.get("tools", {})
    if not isinstance(tools, dict):
        return
    micom_enabled = bool(community.get("enabled", True)) and bool(
        tools.get("micom", True)
    )
    if not micom_enabled:
        return
    binning = config.get("binning", {})
    if not isinstance(binning, dict):
        return
    coverm = binning.get("coverm", {})
    if not isinstance(coverm, dict):
        return
    domains = coverm.get("domains", ["prok", "euk"])
    methods = coverm.get(
        "methods",
        ["relative_abundance", "mean", "covered_fraction", "length", "count"],
    )
    if isinstance(domains, list) and "prok" not in domains:
        raise WorkflowError(
            "Integrated metagenome MICOM requires 'prok' in " "binning.coverm.domains."
        )
    if isinstance(methods, list) and "relative_abundance" not in methods:
        raise WorkflowError(
            "Integrated metagenome MICOM requires 'relative_abundance' in "
            "binning.coverm.methods."
        )


validate_integrated_micom_coverm_config()
validate(config, schema="../schemas/config.schema.yaml")

MAIN_RESOURCES_ROOT = config.get("resources_root", "")


def inherit_main_resources(module_config):
    """Forward a nonblank shared root only when the module root is blank."""
    forwarded = dict(module_config)
    if (
        str(MAIN_RESOURCES_ROOT).strip()
        and not str(forwarded.get("resources_root", "")).strip()
    ):
        forwarded["resources_root"] = MAIN_RESOURCES_ROOT
    return forwarded


if "modules_root" in config:
    config["modules_root"] = expand_configured_path(config["modules_root"])
for section, key in [
    ("porechop_abi", "ab_initio_config"),
    ("porechop_abi", "custom_adapters"),
]:
    if key in config["preprocessing"].get(section, {}):
        config["preprocessing"][section][key] = expand_configured_path(
            config["preprocessing"][section][key]
        )

MODULES_ROOT, MODULE_SNAKEFILES = _resolve_modules_root()

preprocessing_config = config["preprocessing"]
module_preprocessing_config = MODULES_ROOT / "MeeW-Preprocessing" / "config"
if not preprocessing_config.get("porechop_abi", {}).get("ab_initio_config"):
    preprocessing_config.setdefault("porechop_abi", {})["ab_initio_config"] = str(
        module_preprocessing_config / "porechop_abi.config"
    )
