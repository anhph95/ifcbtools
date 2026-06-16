"""Cast aggregation and bottle-merge helpers."""

from __future__ import annotations

import logging
import os
from typing import Sequence

import numpy as np
import pandas as pd

from .constants import DEFAULT_BOTTLE_URL_TEMPLATE
from .stations import StationLocator

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


def aggregate_cast_data(
    df: pd.DataFrame,
    data_cols: Sequence[str],
    station_reference: str | os.PathLike[str] | None = None,
    max_station_distance_km: float | None = 2.0,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    skip_bottle_merge: bool = False,
) -> pd.DataFrame:
    """Aggregate cast replicate samples, assign nearest stations, and merge CTD bottle metadata."""
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

    cruise_frames: list[pd.DataFrame] = []
    bottle_cache: dict[str, pd.DataFrame] = {}
    locator = StationLocator(station_reference=station_reference, max_distance_km=max_station_distance_km)

    for cruise in cast_agg["cruise"].dropna().unique():
        LOGGER.info("Processing cruise: %s", cruise)
        sub = cast_agg.loc[cast_agg["cruise"] == cruise].copy()

        try:
            names, distances = locator.nearest_stations(
                sub,
                lat_col="latitude",
                lon_col="longitude",
                time_col="sample_time",
                max_distance_km=max_station_distance_km,
            )
            sub["nearest_station"] = names
            sub["station_distance"] = distances
        except Exception as exc:
            LOGGER.warning("Station lookup failed for %s: %s", cruise, exc)
            sub["nearest_station"] = pd.NA
            sub["station_distance"] = np.nan

        if skip_bottle_merge:
            cruise_frames.append(sub)
            continue

        try:
            if cruise not in bottle_cache:
                url = bottle_url_template.format(cruise=str(cruise).lower())
                bottle_cache[cruise] = pd.read_csv(url)
                for key in ["cast", "niskin"]:
                    if key in bottle_cache[cruise].columns:
                        bottle_cache[cruise][key] = normalize_integer_key(bottle_cache[cruise][key])

            bottle = bottle_cache[cruise].drop_duplicates(subset=["cast", "niskin"])
            merged = pd.merge(
                sub,
                bottle,
                on=["cast", "niskin"],
                how="left",
                suffixes=("", "_bottle"),
                validate="many_to_one",
            )
        except Exception as exc:
            LOGGER.warning("No CTD bottle data for %s; skipping merge. %s", cruise, exc)
            merged = sub

        cruise_frames.append(merged)

    cast_final = pd.concat(cruise_frames, ignore_index=True).convert_dtypes()
    other_data = other_data.convert_dtypes()
    frames = [frame.dropna(axis=1, how="all") for frame in [cast_final, other_data] if not frame.empty]
    result = pd.concat(frames, ignore_index=True)
    result["sample_time"] = pd.to_datetime(result["sample_time"], errors="coerce")
    return result.sort_values("sample_time").reset_index(drop=True)
