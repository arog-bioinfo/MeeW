import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REAL_MODULES_ROOT = Path(os.environ.get("MEEW_MODULES_ROOT") or ROOT.parent).resolve()


@pytest.fixture
def all_modules_root(tmp_path):
    """Use all five real sibling module sources through a temporary root."""
    modules_root = tmp_path / "scoped_modules"
    modules_root.mkdir()
    for repository in (
        "MeeW-Preprocessing",
        "MeeW-Assembly",
        "MeeW-Binning",
        "MeeW-Annotation",
        "MeeW-Metabolic_Modeling",
    ):
        (modules_root / repository).symlink_to(
            REAL_MODULES_ROOT / repository, target_is_directory=True
        )
    return modules_root
