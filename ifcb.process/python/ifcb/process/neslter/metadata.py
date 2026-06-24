"""Metadata cleaning helpers."""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

DEFAULT_SAMPLE_TYPES = ("cast", "underway", "underway_discrete")


def process_meta(meta: pd.DataFrame, sample_type: Sequence[str] | None = None) -> pd.DataFrame:
    """Clean and enrich IFCB metadata."""
    meta = meta.copy()
    required = ["skip", "sample_type", "depth", "ml_analyzed", "longitude", "latitude", "sample_time"]
    missing = [col for col in required if col not in meta.columns]
    if missing:
        raise ValueError(f"metadata missing required columns: {missing}")

    selected_sample_types = tuple(sample_type) if sample_type else DEFAULT_SAMPLE_TYPES
    meta = meta.loc[(meta["skip"] == 0) & (meta["sample_type"].isin(selected_sample_types))].copy()
    meta.loc[meta["sample_type"] == "underway_discrete", "sample_type"] = "underway"

    meta["depth"] = pd.to_numeric(meta["depth"], errors="coerce")
    mask = meta["sample_type"].isin(["underway", "bucket"])
    meta.loc[mask, "depth"] = meta.loc[mask, "depth"].fillna(0)

    meta["ml_analyzed"] = pd.to_numeric(meta["ml_analyzed"], errors="coerce").replace(0, np.nan)
    meta["longitude"] = pd.to_numeric(meta["longitude"], errors="coerce")
    meta["latitude"] = pd.to_numeric(meta["latitude"], errors="coerce")
    meta = meta.dropna(subset=["sample_time", "longitude", "latitude", "ml_analyzed"])

    meta["sample_time"] = pd.to_datetime(meta["sample_time"], errors="coerce")
    meta = meta.dropna(subset=["sample_time"])
    meta["year"] = meta["sample_time"].dt.year
    meta["month"] = meta["sample_time"].dt.month
    meta["day"] = meta["sample_time"].dt.day
    meta["week"] = meta["sample_time"].dt.isocalendar().week
    meta["doy"] = meta["sample_time"].dt.dayofyear
    meta["season"] = meta["month"].map(
        lambda m: "JFM" if m <= 3 else "AMJ" if m <= 6 else "JAS" if m <= 9 else "OND"
    )

    meta["nearest_station"] = pd.NA
    meta["station_distance"] = np.nan
    return meta.reset_index(drop=True)
