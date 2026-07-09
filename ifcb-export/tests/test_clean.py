"""Tests for the clean dataframe workflow."""

from __future__ import annotations

import unittest

import pandas as pd

from ifcb.clean import filter_and_normalize


class CleanTests(unittest.TestCase):
    """Verify the public clean operation."""

    def make_raw(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "pid": ["underway_1", "discrete_1", "junk_1", "skip_1"],
                "skip": [0, 0, 0, 1],
                "sample_type": ["underway", "underway_discrete", "junk", "underway"],
                "depth": [pd.NA, pd.NA, 0, 0],
                "ml_analyzed": [2.0, 4.0, 1.0, 1.0],
                "longitude": [-70.0, -70.0, -70.0, -70.0],
                "latitude": [40.0, 40.0, 40.0, 40.0],
                "sample_time": [
                    "2020-04-01T00:00:00Z",
                    "2020-07-01T00:00:00Z",
                    "2020-10-01T00:00:00Z",
                    "2020-01-01T00:00:00Z",
                ],
                "taxon_a": [2.0, 8.0, 5.0, 1.0],
            }
        )

    def test_filter_and_normalize_filters_skipped_rows_and_scales_taxa(self) -> None:
        taxonomy = pd.DataFrame({"Annotations": ["taxon_a"]})

        result = filter_and_normalize(self.make_raw(), taxonomy, scaling_factor=1000.0)

        self.assertEqual(result["pid"].tolist(), ["underway_1", "discrete_1", "junk_1"])
        self.assertEqual(result["sample_type"].tolist(), ["underway", "underway_discrete", "junk"])
        self.assertEqual(result["depth"].tolist(), [0, 0, 0])
        self.assertEqual(result["taxon_a"].tolist(), [1000.0, 2000.0, 5000.0])

if __name__ == "__main__":
    unittest.main()
