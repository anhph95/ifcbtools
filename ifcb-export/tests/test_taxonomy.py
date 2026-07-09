"""Tests for consistent taxonomy labels."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from ifcb.taxonomy import import_google_sheet, main, map_taxa_to_label, taxonomy_mapping


class TaxonomyTests(unittest.TestCase):
    """Verify downloaded taxonomy tables use the shared Label column."""

    def test_download_builds_label_from_deepest_taxonomic_level(self) -> None:
        """Create Label from the deepest non-missing level in each row."""
        source = pd.DataFrame(
            {
                "Phylum": ["Bacillariophyta", "Dinoflagellata"],
                "Class": ["Bacillariophyceae", pd.NA],
                "Genus": ["Thalassiosira", pd.NA],
            }
        )

        with patch("ifcb.taxonomy.pd.read_csv", return_value=source):
            taxonomy = import_google_sheet("https://docs.google.com/spreadsheets/d/example/edit?gid=0")

        self.assertEqual(taxonomy["Label"].tolist(), ["Thalassiosira", "Dinoflagellata"])

    def test_taxonomy_mapping_alias_aggregates_columns(self) -> None:
        """Use the public alias to aggregate annotation columns to genus."""
        df = pd.DataFrame(
            {
                "pid": ["sample1", "sample2"],
                "taxon_a": [1, 2],
                "taxon_b": [3, 4],
                "taxon_c": [5, 6],
            }
        )
        taxonomy = pd.DataFrame(
            {
                "Annotations": ["taxon_a", "taxon_b", "taxon_c"],
                "Genus": ["Shared", "Shared", "Other"],
            }
        )

        mapped, taxa = taxonomy_mapping(df, taxonomy)

        self.assertEqual(taxa, ["Other", "Shared"])
        self.assertEqual(mapped["pid"].tolist(), ["sample1", "sample2"])
        self.assertEqual(mapped["Shared"].tolist(), [4, 6])
        self.assertEqual(mapped["Other"].tolist(), [5, 6])

    def test_cli_writes_mapped_csv(self) -> None:
        """Run the CLI entry point against CSV files and write mapped output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            input_path = tmp_path / "input.csv"
            taxonomy_path = tmp_path / "taxonomy.csv"
            output_path = tmp_path / "mapped.csv"
            pd.DataFrame(
                {
                    "pid": ["sample1", "sample2"],
                    "taxon_a": [1, 2],
                    "taxon_b": [3, 4],
                }
            ).to_csv(input_path, index=False)
            pd.DataFrame(
                {
                    "Annotations": ["taxon_a", "taxon_b"],
                    "Genus": ["Shared", "Shared"],
                }
            ).to_csv(taxonomy_path, index=False)

            exit_code = main(
                [
                    "--input",
                    str(input_path),
                    "--taxonomy-file",
                    str(taxonomy_path),
                    "--output",
                    str(output_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            mapped = pd.read_csv(output_path)
            self.assertEqual(mapped["pid"].tolist(), ["sample1", "sample2"])
            self.assertEqual(mapped["Shared"].tolist(), [4, 6])

    def test_map_taxa_to_label_preserves_all_mapped_taxa(self) -> None:
        """Mapping does not hard-code taxon exclusions."""
        df = pd.DataFrame(
            {
                "pid": ["sample1"],
                "nanoplankton_mix": [2],
                "taxon_a": [3],
            }
        )
        taxonomy = pd.DataFrame(
            {
                "Annotations": ["nanoplankton_mix", "taxon_a"],
                "Label": ["Nanoplankton", "Other"],
            }
        )

        mapped, taxa, _ = map_taxa_to_label(df, taxonomy)

        self.assertEqual(taxa, ["Nanoplankton", "Other"])
        self.assertEqual(mapped["Nanoplankton"].tolist(), [2])
        self.assertEqual(mapped["Other"].tolist(), [3])

if __name__ == "__main__":
    unittest.main()
