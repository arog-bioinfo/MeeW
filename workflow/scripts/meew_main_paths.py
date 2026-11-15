"""Portable module-root resolution for the MeeW orchestrator."""

import os
from pathlib import Path

MODULE_REPOSITORIES = {
    "preprocessing": "MeeW-Preprocessing",
    "assembly": "MeeW-Assembly",
    "binning": "MeeW-Binning",
    "annotation": "MeeW-Annotation",
    "metabolic_modeling": "MeeW-Metabolic_Modeling",
}


def expand_path(value):
    """Expand environment variables and ``~`` without rebasing relative paths."""
    if value is None or str(value).strip() == "":
        return ""
    return os.path.expanduser(os.path.expandvars(str(value)))


def workdir_absolute_path(value, cwd=None):
    """Resolve a configured module root relative to the Snakemake workdir."""
    path = Path(expand_path(value))
    if not path.is_absolute():
        path = Path(cwd or Path.cwd()) / path
    return path.resolve()


def resolve_modules_root(configured_root, workflow_basedir, environ=None, cwd=None):
    """Resolve configured, environment, then side-by-side module roots."""
    environ = os.environ if environ is None else environ
    environment_root = environ.get("MEEW_MODULES_ROOT", "")
    if configured_root and str(configured_root).strip():
        root = workdir_absolute_path(configured_root, cwd=cwd)
    elif environment_root.strip():
        root = workdir_absolute_path(environment_root, cwd=cwd)
    else:
        root = Path(workflow_basedir).resolve().parent.parent

    snakefiles = {
        name: root / repository / "workflow" / "Snakefile"
        for name, repository in MODULE_REPOSITORIES.items()
    }
    missing = [str(path) for path in snakefiles.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            f"Module Snakefile(s) missing under {root}: {', '.join(missing)}"
        )
    return root, snakefiles
