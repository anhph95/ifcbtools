"""Tests for IFCB per-liter normalization."""

from __future__ import annotations

import unittest

import pandas as pd

from ifcb.process.neslter.normalize import normalize


class NormalizeTests(unittest.TestCase):
    """Verify one numeric treatment for every selected input file."""

    def test_normalization_preserves_fractional_values(self) -> None:
        """Convert to per-liter values without count-specific rounding."""
        source = pd.DataFrame({"ml_analyzed": [3.0], "taxon": [1.0]})

        result = normalize(source, ["taxon"])

        self.assertAlmostEqual(result.loc[0, "taxon"], 1000 / 3)


if __name__ == "__main__":
    unittest.main()
