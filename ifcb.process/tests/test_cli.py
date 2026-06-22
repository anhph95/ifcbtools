"""Tests for independently selectable IFCB processing operations."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ifcb.process.neslter.cli import parse_args
from ifcb.process.neslter.process import default_output_path, process


class ProcessCliTests(unittest.TestCase):
    """Verify operation selection and enrichment-only processing."""

    def test_cli_requires_at_least_one_operation(self) -> None:
        """Reject commands that do not select any pipeline work."""
        with self.assertRaises(SystemExit):
            parse_args(["input.csv"])

    def test_cli_operations_are_independent(self) -> None:
        """Accept enrichment-only operation selections without --clean."""
        args = parse_args(["data/example.csv", "--add-station", "--merge-nutrient"])

        self.assertFalse(args.clean)
        self.assertTrue(args.add_station)
        self.assertFalse(args.merge_bottle)
        self.assertTrue(args.merge_nutrient)

    def test_clean_requires_one_data_type(self) -> None:
        """Require an explicit normalization type only for cleaning."""
        with self.assertRaises(SystemExit):
            parse_args(["input.csv", "--clean"])

        args = parse_args(["input.csv", "--clean", "--data-type", "count"])
        self.assertEqual(args.data_type, "count")

    def test_default_output_uses_operation_suffixes(self) -> None:
        """Build one cumulative output name in pipeline order."""
        output = default_output_path(
            "data/sample.csv",
            clean=True,
            add_station=True,
            merge_nutrient=True,
        )

        self.assertEqual(output, Path("data/sample_clean_station_nutrient.csv"))

    def test_enrichment_only_reads_existing_output(self) -> None:
        """Apply a selected enrichment operation to exactly one input file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            input_path = output_dir / "sample_clean.csv"
            output_path = output_dir / "sample_clean_nutrient.csv"
            input_path.write_text("pid,value\nsample_1,2\n", encoding="utf-8")

            def fake_merge_nutrients(df, nutrient_source):
                enriched = df.copy()
                enriched["nutrient"] = 3.5
                return enriched

            with patch("ifcb.process.neslter.process.merge_nutrients", side_effect=fake_merge_nutrients):
                outputs = process(
                    input_file=input_path,
                    clean=False,
                    add_station=False,
                    merge_bottle_data=False,
                    merge_nutrient=True,
                )

            self.assertEqual(outputs, output_path)
            self.assertIn("nutrient", output_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
