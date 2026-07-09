"""Tests for the main IFCB process workflow."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from ifcb.process import parse_args, process


class ProcessTests(unittest.TestCase):
    """Verify the public process operation."""

    def test_parse_args_accepts_positional_input(self) -> None:
        args = parse_args(["sample.csv", "--all"])

        self.assertEqual(args.input_file, "sample.csv")
        self.assertTrue(args.clean)
        self.assertTrue(args.station)
        self.assertTrue(args.bottle)
        self.assertTrue(args.nutrient)

    def test_parse_args_accepts_input_flag(self) -> None:
        args = parse_args(["--input", "sample.csv", "--nutrient"])

        self.assertEqual(args.input_file, "sample.csv")
        self.assertTrue(args.nutrient)

    def test_parse_args_rejects_conflicting_inputs(self) -> None:
        with self.assertRaises(SystemExit):
            parse_args(["one.csv", "--input", "two.csv", "--nutrient"])

    def test_enrichment_only_writes_final_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            input_path = output_dir / "sample_clean.csv"
            output_path = output_dir / "sample_clean_nutrient.csv"
            input_path.write_text("pid,value\nsample_1,2\n", encoding="utf-8")

            def fake_nutrient(df, nutrient_source):
                enriched = df.copy()
                enriched["nutrient"] = 3.5
                return enriched

            with patch("ifcb.add.nutrient", side_effect=fake_nutrient):
                output = process(input_file=input_path, nutrient=True)

            self.assertEqual(output, output_path)
            self.assertIn("nutrient", output_path.read_text(encoding="utf-8"))

    def test_clean_and_station_write_only_final_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            input_path = output_dir / "ifcb_count.csv"
            final_path = output_dir / "ifcb_count_clean_station.csv"
            pd.DataFrame(
                {
                    "pid": ["sample_1"],
                    "skip": [0],
                    "sample_type": ["underway"],
                    "depth": [0],
                    "ml_analyzed": [2.0],
                    "longitude": [-70.0],
                    "latitude": [40.0],
                    "sample_time": ["2020-01-01T00:00:00Z"],
                    "taxon_a": [2.0],
                }
            ).to_csv(input_path, index=False)
            pd.DataFrame({"Annotations": ["taxon_a"]}).to_csv(output_dir / "ifcb_taxonomy.csv", index=False)

            def fake_station(df, **kwargs):
                out = df.copy()
                out["nearest_station"] = "L1"
                return out

            with patch("ifcb.add.nearest_station", side_effect=fake_station):
                output = process(input_path, clean=True, station=True)

            self.assertEqual(output, final_path)
            self.assertTrue(final_path.exists())
            self.assertFalse((output_dir / "ifcb_count_clean.csv").exists())
            result = pd.read_csv(final_path)
            self.assertEqual(result.loc[0, "taxon_a"], 1000.0)


if __name__ == "__main__":
    unittest.main()
