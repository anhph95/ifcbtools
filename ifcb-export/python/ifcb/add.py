"""Reusable IFCB post-clean processing operations."""

from __future__ import annotations

import logging
import os
from typing import Sequence

import numpy as np
import pandas as pd

LOGGER = logging.getLogger("ifcb")


def normalize_integer_key(series: pd.Series) -> pd.Series:
    """Extract integer-like values from a merge-key series as nullable Int64."""
    return (
        series.astype(str)
        .str.extract(r"(\d+)", expand=False)
        .pipe(pd.to_numeric, errors="coerce")
        .astype("Int64")
    )


def nearest_station(
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


def bottle(
    df: pd.DataFrame,
    bottle_url_template: str = "https://nes-lter-api.whoi.edu/api/ctd/bottles/{cruise}.csv",
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


DEFAULT_NUTRIENT_COLS = (
    "nitrate_nitrite",
    "ammonium",
    "phosphate",
    "silicate",
    "flag_nitrate_nitrite",
    "flag_ammonium",
    "flag_phosphate",
    "flag_silicate",
)


def fill_cast_from_udw_from_bottles(
    df: pd.DataFrame,
    bottle_url_template: str = "https://nes-lter-api.whoi.edu/api/ctd/bottles/{cruise}.csv",
    sample_type_col: str = "sample_type",
    cruise_col: str = "cruise",
    cast_col: str = "cast",
    depth_col_bottle: str = "depsm",
    niskin_col: str = "niskin",
    target_sample_type: str = "cast_from_udw",
    max_depth: float = 10.0,
    use_median: bool = True,
) -> pd.DataFrame:
    """Fill CTD bottle metadata for synthetic underway-derived cast rows."""
    out = df.copy()
    if niskin_col in out.columns:
        out[niskin_col] = out[niskin_col].astype("object")
    out["_cruise_key"] = out[cruise_col].astype(str).str.strip()
    out["_cast_key"] = pd.to_numeric(out[cast_col], errors="coerce").astype("Int64")
    out["_sample_type_key"] = out[sample_type_col].astype(str).str.strip().str.lower()
    target_mask = out["_sample_type_key"].eq(target_sample_type.lower()) & out["_cast_key"].notna()
    LOGGER.info("Bottle-filling %s %s rows", int(target_mask.sum()), target_sample_type)

    for cruise in out.loc[target_mask, "_cruise_key"].dropna().unique():
        url = bottle_url_template.format(cruise=str(cruise).lower())
        try:
            bottle = pd.read_csv(url, low_memory=False)
        except Exception as exc:
            LOGGER.warning("Skipping bottle fill for %s: %s", cruise, exc)
            continue

        bottle.columns = bottle.columns.str.strip()
        if {cast_col, depth_col_bottle} - set(bottle.columns):
            LOGGER.warning("Skipping bottle fill for %s: required bottle columns are missing.", cruise)
            continue

        bottle = bottle.copy()
        bottle["_cast_key"] = pd.to_numeric(bottle[cast_col], errors="coerce").astype("Int64")
        bottle[depth_col_bottle] = pd.to_numeric(bottle[depth_col_bottle], errors="coerce")
        bottle = bottle.loc[bottle["_cast_key"].notna() & bottle[depth_col_bottle].notna()]
        bottle = bottle.loc[bottle[depth_col_bottle] <= max_depth].copy()
        if bottle.empty:
            continue

        if use_median:
            numeric_cols = [col for col in bottle.select_dtypes(include=[np.number]).columns if col != cast_col]
            bottle_pick = bottle.groupby("_cast_key", as_index=False)[numeric_cols].median()
        else:
            bottle_pick = bottle.sort_values(["_cast_key", depth_col_bottle]).groupby("_cast_key", as_index=False).first()

        if niskin_col in bottle.columns:
            bottle[niskin_col] = bottle[niskin_col].astype("string").str.strip()
            niskin_used = (
                bottle.groupby("_cast_key")[niskin_col]
                .apply(lambda s: ",".join(x for x in s.dropna().unique() if x != ""))
                .reset_index(name="_niskin_used")
            )
            bottle_pick = bottle_pick.merge(niskin_used, on="_cast_key", how="left")

        row_mask = out["_cruise_key"].eq(cruise) & target_mask
        bottle_lookup = bottle_pick.set_index("_cast_key")
        shared_cols = [
            col
            for col in bottle_pick.columns
            if col in out.columns and col not in {cast_col, "_cast_key", "_niskin_used", niskin_col}
        ]

        for idx in out.index[row_mask]:
            cast_key = out.at[idx, "_cast_key"]
            if pd.isna(cast_key) or cast_key not in bottle_lookup.index:
                continue
            src = bottle_lookup.loc[cast_key]
            if isinstance(src, pd.DataFrame):
                src = src.iloc[0]
            for col in shared_cols:
                if pd.isna(out.at[idx, col]) and pd.notna(src[col]):
                    out.at[idx, col] = src[col]
            if niskin_col in out.columns and "_niskin_used" in src.index:
                current = out.at[idx, niskin_col]
                niskin_value = src["_niskin_used"]
                if pd.isna(current) and pd.notna(niskin_value) and str(niskin_value).strip():
                    out.at[idx, niskin_col] = niskin_value
        LOGGER.info("Filled bottle fields for %s %s rows in cruise %s", int(row_mask.sum()), target_sample_type, cruise)

    return out.drop(columns=["_cruise_key", "_cast_key", "_sample_type_key"], errors="ignore")


def nutrient(
    df: pd.DataFrame,
    nutrient_source: str | os.PathLike[str] = "https://nes-lter-api.whoi.edu/api/nut/all.csv",
    nutrient_cols: Sequence[str] = DEFAULT_NUTRIENT_COLS,
    cruise_col: str = "cruise",
    cast_col: str = "cast",
    niskin_col: str = "niskin",
) -> pd.DataFrame:
    """Attach median nutrient values by cruise, cast, and candidate niskins."""
    nut = pd.read_csv(nutrient_source, low_memory=False)
    available_cols = [col for col in nutrient_cols if col in nut.columns]
    if not available_cols:
        LOGGER.warning("No expected nutrient columns found in %s.", nutrient_source)
        return df.copy()

    nut_avg = nut.groupby([cruise_col, niskin_col, cast_col], as_index=False)[available_cols].mean()
    LOGGER.info("Loaded nutrient rows: %s; available nutrient columns: %s", len(nut), len(available_cols))
    nut_avg[cast_col] = pd.to_numeric(nut_avg[cast_col], errors="coerce")
    nut_avg[niskin_col] = pd.to_numeric(nut_avg[niskin_col], errors="coerce")

    out = df.copy()
    out["_row_id"] = np.arange(len(out))
    out["_cast_numeric"] = pd.to_numeric(out[cast_col], errors="coerce")
    out["_niskin_candidates"] = out[niskin_col].astype("string").str.findall(r"\d+")
    exploded = out.explode("_niskin_candidates").copy()
    exploded["_niskin_numeric"] = pd.to_numeric(exploded["_niskin_candidates"], errors="coerce")

    merged = exploded.merge(
        nut_avg,
        left_on=[cruise_col, "_niskin_numeric", "_cast_numeric"],
        right_on=[cruise_col, niskin_col, cast_col],
        how="left",
        suffixes=("", "_nut"),
    )
    source_cols = {col: f"{col}_nut" if f"{col}_nut" in merged.columns else col for col in available_cols}
    matched = merged.loc[merged[list(source_cols.values())].notna().any(axis=1)].copy()
    if matched.empty:
        LOGGER.info("No nutrient matches found")
        return out.drop(columns=["_row_id", "_cast_numeric", "_niskin_candidates"], errors="ignore")

    # For filled rows, niskin may be "21,22,23,24"; median all matching bottles.
    row_medians = matched.groupby("_row_id", as_index=False)[list(source_cols.values())].median()
    LOGGER.info("Merged nutrients for %s IFCB rows using candidate niskins", len(row_medians))
    for col in available_cols:
        source_col = source_cols[col]
        if source_col in row_medians.columns:
            out = out.merge(
                row_medians[["_row_id", source_col]].rename(columns={source_col: f"_{col}_median"}),
                on="_row_id",
                how="left",
            )
            out[col] = out[f"_{col}_median"].combine_first(out[col] if col in out.columns else pd.Series(pd.NA, index=out.index))
            out = out.drop(columns=[f"_{col}_median"])

    return out.drop(columns=["_row_id", "_cast_numeric", "_niskin_candidates"], errors="ignore")
