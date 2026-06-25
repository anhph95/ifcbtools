"""Tests for independently selectable IFCB processing operations."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ifcb.process.neslter.cli import parse_args
from ifcb.process.neslter.process import default_output_path, normalization_settings, process, resolve_data_type


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

    def test_clean_accepts_explicit_input_without_data_type(self) -> None:
        """Use the selected input directly without a redundant type flag."""
        args = parse_args(["input.csv", "--clean"])

        self.assertTrue(args.clean)
        self.assertIsNone(args.data_type)

    def test_clean_accepts_explicit_carbon_data_type(self) -> None:
        """Allow carbon input to select carbon-specific normalization."""
        args = parse_args(["input.csv", "--clean", "--data-type", "carbon"])

        self.assertTrue(args.clean)
        self.assertEqual(args.data_type, "carbon")

    def test_all_enables_every_operation(self) -> None:
        """Expand --all into the complete processing pipeline."""
        args = parse_args(["input.csv", "--all"])

        self.assertTrue(args.all)
        self.assertTrue(args.clean)
        self.assertTrue(args.add_station)
        self.assertTrue(args.merge_bottle)
        self.assertTrue(args.merge_nutrient)

    def test_default_output_uses_operation_suffixes(self) -> None:
        """Build one cumulative output name in pipeline order."""
        output = default_output_path(
            "data/sample.csv",
            clean=True,
            add_station=True,
            merge_nutrient=True,
        )

        self.assertEqual(output, Path("data/sample_clean_station_nutrient.csv"))

    def test_count_data_type_uses_cells_per_liter_scale(self) -> None:
        """Use count normalization for cell abundance products."""
        self.assertEqual(normalization_settings("count"), (1000.0, "cells L-1"))

    def test_carbon_data_type_uses_micrograms_per_liter_scale(self) -> None:
        """Use carbon normalization for biomass products."""
        self.assertEqual(normalization_settings("carbon"), (0.001, "ug C L-1"))

    def test_unknown_data_type_fails(self) -> None:
        """Reject product names without defined unit conversions."""
        with self.assertRaises(ValueError):
            normalization_settings("biovolume")

    def test_count_filename_resolves_count_data_type(self) -> None:
        """Infer count normalization from standard count filenames."""
        self.assertEqual(resolve_data_type("data/ifcb_count.csv"), "count")

    def test_carbon_filename_resolves_carbon_data_type(self) -> None:
        """Infer carbon normalization from standard carbon filenames."""
        self.assertEqual(resolve_data_type("data/ifcb_carbon.csv"), "carbon")

    def test_explicit_data_type_overrides_filename(self) -> None:
        """Honor explicit data type when provided."""
        self.assertEqual(resolve_data_type("data/ifcb_count.csv", data_type="carbon"), "carbon")

    def test_ambiguous_filename_requires_data_type(self) -> None:
        """Ask users to rerun with an explicit data type for ambiguous names."""
        with self.assertRaisesRegex(ValueError, "--data-type count or --data-type carbon"):
            resolve_data_type("data/ifcb.csv")

    def test_all_uses_expanded_operation_suffixes(self) -> None:
        """Name --all output from the operations it expands to."""
        args = parse_args(["data/sample.csv", "--all"])
        output = default_output_path(
            "data/sample.csv",
            clean=args.clean,
            add_station=args.add_station,
            merge_bottle_data=args.merge_bottle,
            merge_nutrient=args.merge_nutrient,
        )

        self.assertEqual(output, Path("data/sample_clean_station_bottle_nutrient.csv"))

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
