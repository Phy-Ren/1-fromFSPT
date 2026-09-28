"""Check CLI meanings, JSON output, and the deprecated comparison flag."""

import json
from pathlib import Path
import subprocess
import sys
import unittest


SCRIPT = Path(__file__).with_name("classify.py")


def invoke(*args):
    return subprocess.run(
        [sys.executable, "-B", str(SCRIPT), *args],
        capture_output=True, text=True, check=False,
    )


class ClassificationCliTests(unittest.TestCase):
    def successful_result(self, *args):
        process = invoke(*args)
        self.assertEqual(process.returncode, 0, process.stderr)
        return json.loads(process.stdout), process.stderr

    def test_relative_classification_retains_universal_group(self):
        result, error = self.successful_result("--orders", "2", "--twist", "1")
        self.assertEqual(error, "")
        self.assertEqual(result["convention"], "relative_oneform")
        self.assertEqual(result["invariant_factors"], [16])
        self.assertEqual(result["order"], 16)
        self.assertFalse(result["universal_subgroup"]["quotiented_in_this_calculation"])

    def test_absolute_bordism_is_labeled_as_mathematical_comparison(self):
        result, error = self.successful_result(
            "--orders", "2", "--twist", "1", "--absolute-bordism"
        )
        self.assertEqual(error, "")
        self.assertEqual(result["convention"], "absolute_bordism_torsion_characters")
        self.assertEqual(result["invariant_factors"], [])
        self.assertEqual(result["order"], 1)
        self.assertIn("Mathematical comparison only", result["mathematical_scope"]["interpretation"])
        self.assertIn("has not been established", result["mathematical_scope"]["interpretation"])

    def test_mixed_orders_preserve_the_group_calculation(self):
        relative, _ = self.successful_result("--orders", "4,6", "--twist", "1,0")
        absolute, _ = self.successful_result(
            "--orders", "4,6", "--twist", "1,0", "--absolute-bordism"
        )
        self.assertEqual(relative["invariant_factors"], [4, 48])
        self.assertEqual(relative["order"], 192)
        self.assertEqual(absolute["invariant_factors"], [12])
        self.assertEqual(absolute["order"], 12)

    def test_zero_twist_agrees_but_calculation_labels_differ(self):
        relative, _ = self.successful_result("--orders", "3,9", "--twist", "0,0")
        absolute, _ = self.successful_result(
            "--orders", "3,9", "--twist", "0,0", "--absolute-bordism"
        )
        self.assertEqual(relative["invariant_factors"], [3, 3, 9])
        self.assertEqual(relative["invariant_factors"], absolute["invariant_factors"])
        self.assertNotEqual(relative["convention"], absolute["convention"])
        self.assertEqual(absolute["universal_subgroup"]["order_in_relative_theory"], 1)
        self.assertFalse(absolute["universal_subgroup"]["quotiented_in_this_calculation"])

    def test_legacy_alias_warns_without_corrupting_json(self):
        legacy, error = self.successful_result(
            "--orders", "8", "--twist", "1", "--full"
        )
        absolute, _ = self.successful_result(
            "--orders", "8", "--twist", "1", "--absolute-bordism"
        )
        self.assertEqual(legacy, absolute)
        self.assertEqual(absolute["invariant_factors"], [2])
        self.assertIn("deprecated", error)
        self.assertIn("not an established", error)
        self.assertNotIn("Traceback", error)

    def test_conflicting_flags_are_rejected(self):
        process = invoke("--orders", "2", "--twist", "1", "--full", "--absolute-bordism")
        self.assertEqual(process.returncode, 2)
        self.assertEqual(process.stdout, "")
        self.assertIn("not allowed with argument", process.stderr)

    def test_invalid_characters_and_lists_are_rejected(self):
        invalid = [
            ("--orders", "3", "--twist", "1"),
            ("--orders", "4,6", "--twist", "1"),
            ("--orders", "0", "--twist", "0"),
            ("--orders", "2", "--twist", "2"),
            ("--orders", "4,,6", "--twist", "1,0"),
            ("--orders", "4.0", "--twist", "0"),
        ]
        for args in invalid:
            with self.subTest(args=args):
                process = invoke(*args)
                self.assertEqual(process.returncode, 2)
                self.assertEqual(process.stdout, "")
                self.assertIn("error:", process.stderr)
                self.assertNotIn("Traceback", process.stderr)


if __name__ == "__main__":
    unittest.main()
