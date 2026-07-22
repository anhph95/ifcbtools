"""Filter and normalize MATLAB-exported IFCB dataframes."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

LOGGER = logging.getLogger("ifcb")


def filter_and_normalize(
    raw: pd.DataFrame,
    taxonomy: pd.DataFrame,
    scaling_factor: float = 1.0,
) -> pd.DataFrame:
    """Filter QC rows, aggregate casts, select taxon columns, and scale values."""
    df = raw.copy()
    required = ["skip", "sample_type", "depth", "ml_analyzed", "longitude", "latitude", "sample_time"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"input data missing required metadata columns: {missing}")

    df = df.loc[df["skip"] == 0].copy()
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    df.loc[df["sample_type"].isin(["underway", "underway_discrete", "bucket"]), "depth"] = df.loc[
        df["sample_type"].isin(["underway", "underway_discrete", "bucket"]), "depth"
    ].fillna(0)
    df["ml_analyzed"] = pd.to_numeric(df["ml_analyzed"], errors="coerce").replace(0, np.nan)
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")
    df["latitude"] = pd.to_numeric(df["latitude"], errors="coerce")
    df = df.dropna(subset=["sample_time", "longitude", "latitude", "ml_analyzed"])
    df["sample_time"] = pd.to_datetime(df["sample_time"], errors="coerce", format="mixed")
    df = df.dropna(subset=["sample_time"])
    df["year"] = df["sample_time"].dt.year
    df["month"] = df["sample_time"].dt.month
    df["day"] = df["sample_time"].dt.day
    df["week"] = df["sample_time"].dt.isocalendar().week
    df["doy"] = df["sample_time"].dt.dayofyear
    df["season"] = df["month"].map(lambda m: "JFM" if m <= 3 else "AMJ" if m <= 6 else "JAS" if m <= 9 else "OND")
    df["nearest_station"] = pd.NA
    df["station_distance"] = np.nan
    df = df.reset_index(drop=True)
    LOGGER.info("Prepared input rows: %s", len(df))

    if "Annotations" not in taxonomy.columns:
        raise ValueError("taxonomy must contain an 'Annotations' column.")
    data_cols = taxonomy["Annotations"].dropna().astype(str).tolist()
    LOGGER.info("Loaded taxonomy annotations: %s", len(data_cols))

    raw_data_cols = [col for col in df.columns if col in data_cols]
    LOGGER.info("Detected %s taxon columns", len(raw_data_cols))

    if "cast" in set(df["sample_type"].dropna().unique()):
        cast_data = df.loc[df["sample_type"] == "cast"].copy()
        other_data = df.loc[df["sample_type"] != "cast"].copy()
        LOGGER.info(
            "Aggregating cast rows by cruise/cast/depth: %s cast rows, %s non-cast rows",
            len(cast_data),
            len(other_data),
        )

        def first_unique_or_join(values: pd.Series) -> object:
            vals = values.dropna().drop_duplicates()
            if vals.empty:
                return pd.NA
            if len(vals) == 1:
                return vals.iloc[0]
            return ", ".join(map(str, vals))

        group_cols = ["cruise", "cast", "depth"]
        sum_cols = [col for col in cast_data.columns if col in raw_data_cols or col in ["ml_analyzed", "n_images"]]
        time_cols = ["sample_time"] if "sample_time" in cast_data.columns else []
        metadata_cols = [col for col in cast_data.columns if col not in set(group_cols + sum_cols + time_cols)]
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
        for key in ["cast", "niskin"]:
            if key in cast_agg.columns:
                cast_agg[key] = pd.to_numeric(cast_agg[key], errors="coerce").astype("Int64")

        frames = [
            frame.dropna(axis=1, how="all")
            for frame in [cast_agg.convert_dtypes(), other_data.convert_dtypes()]
            if not frame.empty
        ]
        df = pd.concat(frames, ignore_index=True)
        df["sample_time"] = pd.to_datetime(df["sample_time"], errors="coerce", format="mixed")
        df = df.sort_values("sample_time").reset_index(drop=True)
        LOGGER.info("Cast aggregation completed: %s rows", len(df))
    else:
        LOGGER.info("No cast data found; skipping cast aggregation.")

    if scaling_factor != 1.0:
        LOGGER.info("Normalizing %s taxon columns with scaling factor %s", len(raw_data_cols), scaling_factor)
    out = df.copy()
    for col in raw_data_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce") / out["ml_analyzed"] * scaling_factor
    return out
