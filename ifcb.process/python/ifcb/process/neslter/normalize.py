"""Data normalization helpers."""

from __future__ import annotations

from typing import Sequence

import pandas as pd


def normalize(
    df: pd.DataFrame,
    data_cols: Sequence[str],
    scaling_factor: float = 1000.0,
) -> pd.DataFrame:
    """Normalize annotations by analyzed volume and a product-specific scale."""
    cols = [col for col in data_cols if col in df.columns]
    df = df.copy()
    df[cols] = df[cols].div(df["ml_analyzed"], axis=0).mul(scaling_factor)
    return df.reset_index(drop=True)
