import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import sys

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from meew_main_paths import MODULE_REPOSITORIES, expand_path, resolve_modules_root


class ModuleRootResolutionTests(unittest.TestCase):
    def create_layout(self, root):
        for repository in MODULE_REPOSITORIES.values():
            snakefile = root / repository / "workflow" / "Snakefile"
            snakefile.parent.mkdir(parents=True)
            snakefile.touch()

    def test_precedence_and_side_by_side_fallback(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            configured = temporary / "configured modules"
            environment = temporary / "environment modules"
            fallback = temporary / "fallback"
            for root in [configured, environment, fallback]:
                self.create_layout(root)

            workflow_basedir = fallback / "MeeW" / "workflow"
            workflow_basedir.mkdir(parents=True)
            env = {"MEEW_MODULES_ROOT": str(environment)}

            root, _ = resolve_modules_root(configured, workflow_basedir, env)
            self.assertEqual(root, configured.resolve())
            root, _ = resolve_modules_root("", workflow_basedir, env)
            self.assertEqual(root, environment.resolve())
            root, _ = resolve_modules_root("", workflow_basedir, {})
            self.assertEqual(root, fallback.resolve())

    def test_relative_root_uses_workdir_and_expansion_preserves_relative_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            relative_root = temporary / "modules with spaces"
            self.create_layout(relative_root)
            root, _ = resolve_modules_root(
                "modules with spaces",
                temporary / "unused" / "workflow",
                {},
                cwd=temporary,
            )
            self.assertEqual(root, relative_root.resolve())
            self.assertEqual(expand_path("relative/input.tsv"), "relative/input.tsv")

    def test_tilde_and_environment_expansion(self):
        with patch.dict(
            os.environ,
            {"HOME": "/tmp/test home", "MEEW_TEST_PATH": "resource folder"},
            clear=True,
        ):
            self.assertEqual(
                expand_path("~/$MEEW_TEST_PATH"), "/tmp/test home/resource folder"
            )


if __name__ == "__main__":
    unittest.main()
