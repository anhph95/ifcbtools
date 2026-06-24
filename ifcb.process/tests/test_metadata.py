"""Tests for IFCB metadata cleaning defaults."""

from __future__ import annotations

import unittest

import pandas as pd

from ifcb.process.neslter.metadata import process_meta


class MetadataTests(unittest.TestCase):
    """Verify the always-on clean metadata filters."""

    def make_meta(self) -> pd.DataFrame:
        """Build minimal metadata rows that exercise sample-type filtering."""
        return pd.DataFrame(
            {
                "pid": ["cast_1", "underway_1", "discrete_1", "junk_1", "skip_1"],
                "skip": [0, 0, 0, 0, 1],
                "sample_type": ["cast", "underway", "underway_discrete", "junk", "cast"],
                "depth": [5, pd.NA, pd.NA, 0, 5],
                "ml_analyzed": [1, 1, 1, 1, 1],
                "longitude": [-70, -70, -70, -70, -70],
                "latitude": [40, 40, 40, 40, 40],
                "sample_time": [
                    "2020-01-01 00:00:00+00:00",
                    "2020-04-01 00:00:00+00:00",
                    "2020-07-01 00:00:00+00:00",
                    "2020-10-01 00:00:00+00:00",
                    "2020-01-02 00:00:00+00:00",
                ],
            }
        )

    def test_default_keeps_core_sample_types(self) -> None:
        """Default cleaning keeps cast and underway records only."""
        result = process_meta(self.make_meta())

        self.assertEqual(result["pid"].tolist(), ["cast_1", "underway_1", "discrete_1"])
        self.assertEqual(result["sample_type"].tolist(), ["cast", "underway", "underway"])
        self.assertEqual(result.loc[result["pid"] == "discrete_1", "depth"].iloc[0], 0)

    def test_explicit_sample_type_still_limits_rows(self) -> None:
        """User-selected sample types narrow the always-on skip filter."""
        result = process_meta(self.make_meta(), sample_type=["cast"])

        self.assertEqual(result["pid"].tolist(), ["cast_1"])
        self.assertEqual(result["sample_type"].tolist(), ["cast"])


if __name__ == "__main__":
    unittest.main()
