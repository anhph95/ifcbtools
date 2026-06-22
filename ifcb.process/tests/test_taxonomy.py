"""Tests for consistent taxonomy labels."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from ifcb.process.neslter.fill import make_filled_dataset
from ifcb.process.neslter.taxonomy import import_google_sheet


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

        with patch("ifcb.process.neslter.taxonomy.pd.read_csv", return_value=source):
            taxonomy = import_google_sheet("https://docs.google.com/spreadsheets/d/example/edit?gid=0")

        self.assertEqual(taxonomy["Label"].tolist(), ["Thalassiosira", "Dinoflagellata"])

    def test_fill_does_not_create_taxonomy_output(self) -> None:
        """Keep taxonomy in memory and write only the selected data product."""
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir)
            pd.DataFrame(
                {
                    "Annotations": ["taxon_a"],
                    "Phylum": ["Bacillariophyta"],
                    "Label": ["Bacillariophyta"],
                }
            ).to_csv(data_dir / "ifcb_taxonomy.csv", index=False)
            pd.DataFrame(
                {
                    "sample_time": ["2020-01-01T00:00:00Z"],
                    "taxon_a": [1],
                }
            ).to_csv(data_dir / "ifcb_count_clean.csv", index=False)

            def keep_rows(df, **kwargs):
                out = df.copy()
                out["_fill_created"] = False
                return out

            with (
                patch(
                    "ifcb.process.neslter.fill.map_taxa_to_label",
                    side_effect=lambda df, taxonomy: (df, ["taxon_a"], taxonomy),
                ),
                patch("ifcb.process.neslter.fill.fill_missing_casts_from_underway", side_effect=keep_rows),
            ):
                outputs = make_filled_dataset(data_dir, data_types=["count"])

            self.assertEqual(outputs, [data_dir / "ifcb_count_fill.csv"])
            self.assertEqual(list(data_dir.glob("ifcb_taxonomy*.csv")), [data_dir / "ifcb_taxonomy.csv"])


if __name__ == "__main__":
    unittest.main()
