"""Data normalization helpers."""

from __future__ import annotations

from typing import Sequence

import pandas as pd


def normalize(df: pd.DataFrame, data_cols: Sequence[str], data_type: str) -> pd.DataFrame:
    """Normalize annotation columns by ml analyzed and convert to per-liter values."""
    cols = [col for col in data_cols if col in df.columns]
    df = df.copy()
    df[cols] = df[cols].div(df["ml_analyzed"], axis=0).mul(1000)
    if data_type == "count":
        df[cols] = df[cols].round().astype("Int64")
    return df.reset_index(drop=True)
