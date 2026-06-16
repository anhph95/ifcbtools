"""Cast aggregation helpers."""

from __future__ import annotations

import logging
from typing import Sequence

import pandas as pd

LOGGER = logging.getLogger("ifcb.neslter")


def first_unique_or_join(values: pd.Series) -> object:
    """Return one unique non-null value, joined unique values, or NA."""
    vals = values.dropna().drop_duplicates()
    if vals.empty:
        return pd.NA
    if len(vals) == 1:
        return vals.iloc[0]
    return ", ".join(map(str, vals))


def normalize_integer_key(series: pd.Series) -> pd.Series:
    """Extract integer-like values from a merge-key series as nullable Int64."""
    return (
        series.astype(str)
        .str.extract(r"(\d+)", expand=False)
        .pipe(pd.to_numeric, errors="coerce")
        .astype("Int64")
    )


def aggregate_cast_data(df: pd.DataFrame, data_cols: Sequence[str]) -> pd.DataFrame:
    """Aggregate cast replicate samples without station, bottle, or nutrient joins."""
    if "cast" not in set(df["sample_type"].dropna().unique()):
        LOGGER.info("No cast data found; skipping cast aggregation.")
        return df

    cast_data = df[df["sample_type"] == "cast"]
    other_data = df[df["sample_type"] != "cast"]

    agg_dict = {
        col: (
            "sum"
            if col in data_cols or col in ["ml_analyzed", "n_images"]
            else "min"
            if col == "sample_time"
            else first_unique_or_join
        )
        for col in cast_data.columns
        if col not in ["cruise", "cast", "depth"]
    }

    cast_agg = cast_data.groupby(["cruise", "cast", "depth"], as_index=False).agg(agg_dict)
    for key in ["cast", "niskin"]:
        if key in cast_agg.columns:
            cast_agg[key] = pd.to_numeric(cast_agg[key], errors="coerce").astype("Int64")

    cast_final = cast_agg.convert_dtypes()
    other_data = other_data.convert_dtypes()
    frames = [frame.dropna(axis=1, how="all") for frame in [cast_final, other_data] if not frame.empty]
    result = pd.concat(frames, ignore_index=True)
    result["sample_time"] = pd.to_datetime(result["sample_time"], errors="coerce")
    return result.sort_values("sample_time").reset_index(drop=True)
