import ast
import re
import unittest
from pathlib import Path


class ShellPathSafetyTests(unittest.TestCase):
    def test_snakemake_shell_paths_are_quoted_across_all_repositories(self):
        repositories = [
            Path(__file__).resolve().parents[1],
            *[
                Path(__file__).resolve().parents[2] / name
                for name in [
                    "MeeW-Preprocessing",
                    "MeeW-Assembly",
                    "MeeW-Binning",
                    "MeeW-Annotation",
                    "MeeW-Metabolic_Modeling",
                ]
            ],
        ]
        unsafe_placeholder = re.compile(
            r"\{(?:snakemake\.)?(?:input|output)"
            r"(?:\[[^]]+\]|\.[A-Za-z_][A-Za-z0-9_]*)?\}"
        )
        findings = []
        for repository in repositories:
            for path in repository.glob("workflow/**/*.smk"):
                for line_number, line in enumerate(
                    path.read_text(encoding="utf-8").splitlines(), 1
                ):
                    if unsafe_placeholder.search(line):
                        findings.append(f"{path}:{line_number}:{line.strip()}")
            for path in repository.glob("workflow/**/*.py"):
                source = path.read_text(encoding="utf-8")
                tree = ast.parse(source, filename=str(path))
                for node in ast.walk(tree):
                    if not isinstance(node, ast.Call) or not node.args:
                        continue
                    is_shell_call = (
                        isinstance(node.func, ast.Name) and node.func.id == "shell"
                    ) or (
                        isinstance(node.func, ast.Attribute)
                        and node.func.attr == "shell"
                    )
                    if not is_shell_call:
                        continue
                    command = ast.get_source_segment(source, node.args[0]) or ""
                    for offset, line in enumerate(command.splitlines()):
                        if unsafe_placeholder.search(line):
                            findings.append(
                                f"{path}:{node.args[0].lineno + offset}:{line.strip()}"
                            )
        self.assertEqual(
            findings, [], "Unsafe shell path placeholders:\n" + "\n".join(findings)
        )


if __name__ == "__main__":
    unittest.main()
