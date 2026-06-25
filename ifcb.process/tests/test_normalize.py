"""Tests for IFCB per-liter normalization."""

from __future__ import annotations

import unittest

import pandas as pd

from ifcb.process.neslter.normalize import normalize


class NormalizeTests(unittest.TestCase):
    """Verify product-specific per-liter normalization."""

    def test_count_normalization_uses_cells_per_liter_scale(self) -> None:
        """Convert cell counts to cells per liter."""
        source = pd.DataFrame({"ml_analyzed": [3.0], "taxon": [1.0]})

        result = normalize(source, ["taxon"])

        self.assertAlmostEqual(result.loc[0, "taxon"], 1000 / 3)

    def test_carbon_normalization_uses_micrograms_per_liter_scale(self) -> None:
        """Convert pg C product values to ug C per liter."""
        source = pd.DataFrame({"ml_analyzed": [3.0], "taxon": [1.0]})

        result = normalize(source, ["taxon"], scaling_factor=0.001)

        self.assertAlmostEqual(result.loc[0, "taxon"], 0.001 / 3)


if __name__ == "__main__":
    unittest.main()
