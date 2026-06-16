"""Cast aggregation helpers."""

from __future__ import annotations

import logging
from typing import Sequence

import pandas as pd

LOGGER = logging.getLogger("ifcb.process.neslter")


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

    cast_data = df.loc[df["sample_type"] == "cast"].copy()
    other_data = df.loc[df["sample_type"] != "cast"].copy()
    cast_data = cast_data.copy()
    LOGGER.info(
        "Aggregating cast rows by cruise/cast/depth: %s cast rows, %s non-cast rows",
        len(cast_data),
        len(other_data),
    )
    group_cols = ["cruise", "cast", "depth"]
    sum_cols = [col for col in cast_data.columns if col in data_cols or col in ["ml_analyzed", "n_images"]]
    time_cols = ["sample_time"] if "sample_time" in cast_data.columns else []
    metadata_cols = [
        col
        for col in cast_data.columns
        if col not in set(group_cols + sum_cols + time_cols)
    ]

    grouped = cast_data.groupby(group_cols, as_index=False, sort=False, dropna=False)
    pieces = []
    if sum_cols:
        pieces.append(grouped[sum_cols].sum(min_count=1))
    if time_cols:
        pieces.append(grouped[time_cols].min())
    if metadata_cols:
        pieces.append(grouped[metadata_cols].agg(first_unique_or_join))

    cast_agg = pieces[0]
    for piece in pieces[1:]:
        cast_agg = cast_agg.merge(piece, on=group_cols, how="left")
    cast_agg = cast_agg.copy()
    for key in ["cast", "niskin"]:
        if key in cast_agg.columns:
            cast_agg[key] = pd.to_numeric(cast_agg[key], errors="coerce").astype("Int64")

    cast_final = cast_agg.convert_dtypes()
    other_data = other_data.convert_dtypes()
    frames = [frame.dropna(axis=1, how="all") for frame in [cast_final, other_data] if not frame.empty]
    result = pd.concat(frames, ignore_index=True)
    result["sample_time"] = pd.to_datetime(result["sample_time"], errors="coerce")
    result = result.sort_values("sample_time").reset_index(drop=True)
    LOGGER.info("Cast aggregation completed: %s rows", len(result))
    return result
