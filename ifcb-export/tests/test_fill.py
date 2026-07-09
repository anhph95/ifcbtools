"""Tests for the main fill workflow."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import pandas as pd

from ifcb.fill import fill_dataset


class FillTests(unittest.TestCase):
    """Verify the public fill operation."""

    def test_fill_dataset_assigns_station_when_absent(self) -> None:
        taxonomy = pd.DataFrame(
            {
                "Annotations": ["taxon_a"],
                "Phylum": ["Bacillariophyta"],
                "Label": ["Bacillariophyta"],
            }
        )
        df = pd.DataFrame(
            {
                "sample_time": ["2020-01-01T00:00:00Z"],
                "sample_type": ["underway"],
                "latitude": [41.0],
                "longitude": [-70.0],
                "taxon_a": [1],
            }
        )

        def add_station(df, **kwargs):
            out = df.copy()
            out["nearest_station"] = "L1"
            out["station_distance"] = 1.0
            return out

        def keep_rows(df, **kwargs):
            out = df.copy()
            out["_fill_created"] = False
            return out

        with (
            patch("ifcb.fill.load_station_reference", return_value=pd.DataFrame({"station": ["L1", "d4a"]})),
            patch("ifcb.fill.assign_nearest_stations", side_effect=add_station) as assign,
            patch("ifcb.fill.map_taxa_to_label", side_effect=lambda df, taxonomy: (df, ["taxon_a"], taxonomy)),
            patch("ifcb.fill.fill_missing_casts_from_underway", side_effect=keep_rows),
        ):
            output = fill_dataset(df, taxonomy)

        assign.assert_called_once()
        assigned_reference = assign.call_args.kwargs["station_reference"]
        self.assertEqual(assigned_reference["station"].tolist(), ["L1"])
        self.assertEqual(output["nearest_station"].tolist(), ["L1"])


if __name__ == "__main__":
    unittest.main()
