"""Reusable IFCB post-clean processing operations."""

from __future__ import annotations

import logging
import os

import pandas as pd

from .casts import normalize_integer_key
from .constants import DEFAULT_BOTTLE_URL_TEMPLATE
from .fill import fill_cast_from_udw_from_bottles

LOGGER = logging.getLogger("ifcb.process.neslter")


def add_nearest_station(
    df: pd.DataFrame,
    station_reference: pd.DataFrame | str | os.PathLike[str] | None = None,
    max_station_distance_km: float | None = 2.0,
    show_progress: bool = True,
) -> pd.DataFrame:
    """Assign nearest_station and station_distance to every row using StationLocator."""
    from .stations import StationLocator

    LOGGER.info("Assigning nearest stations for %s rows", len(df))
    locator = StationLocator(station_reference=station_reference, max_distance_km=max_station_distance_km)
    out = df.copy()
    out["sample_time"] = pd.to_datetime(out["sample_time"], errors="coerce", utc=True)

    names, distances = locator.nearest_stations(
        out,
        lat_col="latitude",
        lon_col="longitude",
        time_col="sample_time",
        max_distance_km=max_station_distance_km,
        show_progress=show_progress,
    )
    out["nearest_station"] = names
    out["station_distance"] = distances
    LOGGER.info("Assigned nearest stations for %s rows", out["nearest_station"].notna().sum())
    return out


def merge_bottle(
    df: pd.DataFrame,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    cruise_col: str = "cruise",
    cast_col: str = "cast",
    niskin_col: str = "niskin",
) -> pd.DataFrame:
    """Merge CTD bottle data by cruise, cast, and niskin where possible."""
    LOGGER.info("Filling bottle fields for cast_from_udw rows before normal bottle merge")
    df = fill_cast_from_udw_from_bottles(df, bottle_url_template=bottle_url_template)
    out_frames = []
    for cruise, sub in df.groupby(cruise_col, dropna=False, sort=False):
        sub = sub.copy()
        if pd.isna(cruise):
            out_frames.append(sub)
            continue

        try:
            bottle = pd.read_csv(bottle_url_template.format(cruise=str(cruise).lower()), low_memory=False)
        except Exception as exc:
            LOGGER.warning("No CTD bottle data for %s; keeping rows unchanged. %s", cruise, exc)
            out_frames.append(sub)
            continue

        bottle = bottle.copy()
        bottle.columns = bottle.columns.str.strip()

        if cast_col not in bottle.columns or niskin_col not in bottle.columns:
            LOGGER.warning("Bottle data for %s lacks cast/niskin keys; keeping rows unchanged.", cruise)
            out_frames.append(sub)
            continue

        # Use temporary numeric keys so comma-list filled niskins are preserved.
        sub["_cast_key"] = normalize_integer_key(sub[cast_col])
        sub["_niskin_key"] = normalize_integer_key(sub[niskin_col])
        bottle["_cast_key"] = normalize_integer_key(bottle[cast_col])
        bottle["_niskin_key"] = normalize_integer_key(bottle[niskin_col])
        bottle = bottle.drop_duplicates(subset=["_cast_key", "_niskin_key"])
        merged = pd.merge(
            sub,
            bottle,
            on=["_cast_key", "_niskin_key"],
            how="left",
            suffixes=("", "_bottle"),
            validate="many_to_one",
        ).drop(columns=["_cast_key", "_niskin_key"], errors="ignore")
        LOGGER.info("Merged bottle data for cruise %s: %s rows", cruise, len(merged))
        out_frames.append(merged)

    out = pd.concat(out_frames, ignore_index=True).sort_values("sample_time").reset_index(drop=True)
    LOGGER.info("Bottle merge completed: %s rows", len(out))
    return out
