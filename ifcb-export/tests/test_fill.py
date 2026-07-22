"""Tests for the main fill workflow."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import pandas as pd

from ifcb.fill import fill_dataset


class FillTests(unittest.TestCase):
    """Verify the public fill operation."""

    def test_fill_dataset_requires_existing_station_assignment(self) -> None:
        df = pd.DataFrame(
            {
                "sample_time": ["2020-01-01T00:00:00Z"],
                "sample_type": ["underway"],
                "latitude": [41.0],
                "longitude": [-70.0],
                "taxon_a": [1],
            }
        )

        with self.assertRaisesRegex(ValueError, "nearest_station"):
            fill_dataset(df)

    def test_fill_dataset_uses_existing_station_assignment(self) -> None:
        df = pd.DataFrame(
            {
                "sample_time": ["2020-01-01T00:00:00Z"],
                "sample_type": ["underway"],
                "latitude": [41.0],
                "longitude": [-70.0],
                "nearest_station": ["L1"],
                "station_distance": [1.0],
                "taxon_a": [1],
            }
        )

        def keep_rows(df, **kwargs):
            out = df.copy()
            out["_fill_created"] = False
            return out

        with patch("ifcb.fill.fill_missing_casts_from_underway", side_effect=keep_rows):
            output = fill_dataset(df)

        self.assertEqual(output["nearest_station"].tolist(), ["L1"])


if __name__ == "__main__":
    unittest.main()
