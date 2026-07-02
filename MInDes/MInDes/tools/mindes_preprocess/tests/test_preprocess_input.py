from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).parents[1] / "preprocess_input.py"
SPEC = importlib.util.spec_from_file_location("preprocess_input", SCRIPT)
assert SPEC and SPEC.loader
preprocess_input = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(preprocess_input)


class PreprocessInputTests(unittest.TestCase):
    def test_cli_exit_codes_distinguish_no_include_and_input_error(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plain = root / "plain.mindes"
            plain.write_text("A = 1\n", encoding="utf-8")
            malformed = root / "malformed.mindes"
            malformed.write_text("INCLUDE\n", encoding="utf-8")
            output = root / "combined.mindes"

            no_include = subprocess.run(
                [sys.executable, str(SCRIPT), "--input", str(plain), "--output", str(output)],
                check=False,
                capture_output=True,
                text=True,
            )
            input_error = subprocess.run(
                [sys.executable, str(SCRIPT), "--input", str(malformed), "--output", str(output)],
                check=False,
                capture_output=True,
                text=True,
            )

            self.assertEqual(no_include.returncode, 10)
            self.assertEqual(input_error.returncode, 20)
            self.assertIn("INCLUDE requires a file path", input_error.stderr)

    def test_no_include_keeps_original_and_does_not_write_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "plain.mindes"
            output = root / "missing" / "combined_infile.mindes"
            source.write_text("A = 1\n", encoding="utf-8")

            self.assertFalse(preprocess_input.preprocess(source, output))
            self.assertFalse(output.exists())

    def test_nested_include_is_depth_first_and_include_once(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output_dir = root / "case"
            output_dir.mkdir()
            parts = root / "parts"
            parts.mkdir()
            (parts / "leaf.mindes").write_text("LEAF = 1\n", encoding="utf-8")
            (parts / "child.mindes").write_text(
                "CHILD = 1\nINCLUDE leaf.mindes\n", encoding="utf-8"
            )
            source = root / "case.mindes"
            source.write_text(
                "ROOT = before\nINCLUDE parts/child.mindes\n"
                "INCLUDE parts/leaf.mindes\nROOT = after\n",
                encoding="utf-8",
            )
            output = output_dir / "combined_infile.mindes"

            self.assertTrue(preprocess_input.preprocess(source, output))
            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "ROOT = before\nCHILD = 1\nLEAF = 1\nROOT = after\n",
            )

    def test_quoted_unicode_path_bom_and_crlf(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output_dir = root / "入口"
            output_dir.mkdir()
            child = root / "含 空格.mindes"
            child.write_text("值 = 2\r\n", encoding="utf-8-sig", newline="")
            source = root / "入口.mindes"
            source.write_text('  INCLUDE "含 空格.mindes"\r\n', encoding="utf-8-sig", newline="")
            output = output_dir / "combined_infile.mindes"

            self.assertTrue(preprocess_input.preprocess(source, output))
            self.assertEqual(output.read_text(encoding="utf-8"), "值 = 2\n")

    def test_cycle_reports_the_include_chain(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.mindes"
            second = root / "second.mindes"
            first.write_text("INCLUDE second.mindes\n", encoding="utf-8")
            second.write_text("INCLUDE first.mindes\n", encoding="utf-8")

            with self.assertRaisesRegex(preprocess_input.PreprocessError, "circular INCLUDE"):
                preprocess_input.preprocess(first, root / "combined.mindes")

    def test_missing_and_malformed_include_are_errors(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            missing = root / "missing.mindes"
            missing.write_text("INCLUDE absent.mindes\n", encoding="utf-8")
            malformed = root / "malformed.mindes"
            malformed.write_text('INCLUDE "unfinished.mindes\n', encoding="utf-8")

            with self.assertRaises(preprocess_input.PreprocessError):
                preprocess_input.preprocess(missing, root / "missing.out")
            with self.assertRaisesRegex(preprocess_input.PreprocessError, "unterminated"):
                preprocess_input.preprocess(malformed, root / "malformed.out")


if __name__ == "__main__":
    unittest.main()
