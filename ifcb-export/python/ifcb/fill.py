"""Fill missing IFCB cast samples from nearby underway observations."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
import sys
from typing import Sequence

import numpy as np
import pandas as pd

from .add import fill_cast_from_udw_from_bottles, nutrient
from ifcb.logging import log_run_configuration, redact_command_line, setup_logging

LOGGER = logging.getLogger("ifcb")

DEFAULT_TARGET_STATIONS = ("MVCO", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11")


def fill_missing_casts_from_underway(
    df: pd.DataFrame,
    target_stations: Sequence[str] = DEFAULT_TARGET_STATIONS,
    cruise_col: str = "cruise",
    station_col: str = "nearest_station",
    cast_col: str = "cast",
    time_col: str = "sample_time",
    depth_col: str = "depth",
    sample_type_col: str = "sample_type",
    station_distance_col: str = "station_distance",
    cast_label: str = "cast",
    underway_label: str = "underway",
    fill_label: str = "cast_from_udw",
    surface_max_depth: float = 10.0,
    mark_filled_col: str | None = None,
) -> pd.DataFrame:
    """Create synthetic surface casts from underway rows when station casts are missing.

    Ecological goal: preserve a balanced local-community set i across surveys by
    filling absent surface observations with the nearest underway sample from
    the same cruise and station.
    """
    df = df.copy()
    if mark_filled_col is not None:
        df[mark_filled_col] = False
    df[time_col] = pd.to_datetime(df[time_col], errors="coerce", utc=True)
    df[depth_col] = pd.to_numeric(df[depth_col], errors="coerce")
    if station_distance_col in df.columns:
        df[station_distance_col] = pd.to_numeric(df[station_distance_col], errors="coerce")

    # Underway rows are surface samples; blank depth values represent z = 0 m.
    underway_mask = df[sample_type_col].eq(underway_label)
    df.loc[underway_mask, depth_col] = df.loc[underway_mask, depth_col].fillna(0)

    target_stations = tuple(target_stations)
    df = df.loc[df[sample_type_col].isin([cast_label, underway_label])].copy()
    cast_df = df.loc[(df[sample_type_col] == cast_label) & df[station_col].isin(target_stations)].copy()
    udw_df = df.loc[(df[sample_type_col] == underway_label) & df[station_col].isin(target_stations)].copy()
    if udw_df.empty:
        return df.sort_values(time_col).reset_index(drop=True)

    filled = []

    # Rule 1: a cast exists, but its minimum sampled depth is deeper than z = surface_max_depth.
    if not cast_df.empty:
        summary = (
            cast_df.groupby([cruise_col, station_col, cast_col], dropna=False)
            .agg(min_depth=(depth_col, "min"), ref_time=(time_col, "median"))
            .reset_index()
        )
        missing_surface = summary.loc[summary["min_depth"].isna() | (summary["min_depth"] > surface_max_depth)]
        for _, group in missing_surface.iterrows():
            candidates = udw_df.loc[
                (udw_df[cruise_col] == group[cruise_col]) & (udw_df[station_col] == group[station_col])
            ].copy()
            if candidates.empty:
                continue
            candidates["_dt"] = (candidates[time_col] - group["ref_time"]).abs()
            sort_cols = ["_dt"] + ([station_distance_col] if station_distance_col in candidates.columns else [])
            best = candidates.sort_values(sort_cols).iloc[0].copy()
            best[sample_type_col] = fill_label
            best[cast_col] = group[cast_col]
            best[station_col] = group[station_col]
            best[depth_col] = 0
            if mark_filled_col is not None:
                best[mark_filled_col] = True
            filled.append(best.drop(labels=["_dt"], errors="ignore"))

    # Rule 2: a cruise has no cast for a target station; use closest underway for that station.
    cruises = df[cruise_col].dropna().unique()
    existing = cast_df[[cruise_col, station_col]].drop_duplicates()
    expected = pd.MultiIndex.from_product([cruises, target_stations], names=[cruise_col, station_col]).to_frame(
        index=False
    )
    missing_station = expected.merge(existing.assign(has_cast=1), on=[cruise_col, station_col], how="left")
    missing_station = missing_station.loc[missing_station["has_cast"].isna()]
    for _, group in missing_station.iterrows():
        candidates = udw_df.loc[
            (udw_df[cruise_col] == group[cruise_col]) & (udw_df[station_col] == group[station_col])
        ].copy()
        if candidates.empty:
            continue
        sort_cols = ([station_distance_col] if station_distance_col in candidates.columns else []) + [time_col]
        best = candidates.sort_values(sort_cols).iloc[0].copy()
        best[sample_type_col] = fill_label
        best[cast_col] = group[station_col]
        best[station_col] = group[station_col]
        best[depth_col] = 0
        if mark_filled_col is not None:
            best[mark_filled_col] = True
        filled.append(best)

    if filled:
        df = pd.concat([df, pd.DataFrame(filled)], ignore_index=True)
    LOGGER.info("Created %s %s rows from underway data", len(filled), fill_label)
    out = df.sort_values([cruise_col, station_col, cast_col, time_col, depth_col]).reset_index(drop=True)
    if mark_filled_col is None:
        out = out.drop(columns=["_fill_created"], errors="ignore")
    return out


def fill_dataset(
    df: pd.DataFrame,
    target_stations: Sequence[str] = DEFAULT_TARGET_STATIONS,
    bottle_url_template: str = "https://nes-lter-api.whoi.edu/api/ctd/bottles/{cruise}.csv",
    nutrient_source: str | os.PathLike[str] = "https://nes-lter-api.whoi.edu/api/nut/all.csv",
) -> pd.DataFrame:
    """Create missing IFCB cast samples from one dataframe."""
    df = df.copy()
    input_rows = len(df)
    if "nearest_station" not in df.columns:
        raise ValueError("fill input must include nearest_station; run station assignment before fill.")

    df = fill_missing_casts_from_underway(
        df,
        target_stations=target_stations,
        mark_filled_col="_fill_created",
    )
    new_rows = df.loc[df["_fill_created"]].drop(columns=["_fill_created"]).copy()
    existing_rows = df.loc[~df["_fill_created"]].drop(columns=["_fill_created"]).copy()
    if not new_rows.empty:
        new_rows = fill_cast_from_udw_from_bottles(new_rows, bottle_url_template=bottle_url_template)
        new_rows = nutrient(new_rows, nutrient_source=nutrient_source)
    df = pd.concat([existing_rows, new_rows], ignore_index=True)
    out = df.sort_values("sample_time").reset_index(drop=True)
    LOGGER.info("Filled data: %s -> %s rows", input_rows, len(out))
    return out


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create one filled IFCB CSV from one explicitly selected input CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", dest="input_file", required=True, help="Input IFCB CSV file.")
    parser.add_argument("--output", "--output-file", dest="output_file", default=None, help="Output CSV file; defaults to an _fill suffix.")
    parser.add_argument("--bottle-url-template", default="https://nes-lter-api.whoi.edu/api/ctd/bottles/{cruise}.csv")
    parser.add_argument("--nutrient-source", default="https://nes-lter-api.whoi.edu/api/nut/all.csv")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    parser.add_argument("--log-dir", default=None, help="Directory for timestamped workflow log files.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input_file)
    if args.output_file is None:
        output_path = input_path.with_name(f"{input_path.stem}_fill{input_path.suffix}")
    else:
        output_path = Path(args.output_file)
    log_dir = Path(args.log_dir) if args.log_dir is not None else Path.cwd() / "logs"
    setup_logging(log_dir=log_dir, name="ifcb_fill_missing", level=getattr(logging, args.log_level))
    LOGGER.info("Starting IFCB missing-cast fill for %s", input_path)
    log_run_configuration(
        LOGGER,
        {
            "command": redact_command_line(sys.argv),
            "input_file": input_path.resolve(),
            "output_file": output_path.resolve(),
            "bottle_url_template": args.bottle_url_template,
            "nutrient_source": args.nutrient_source,
            "log_level": args.log_level,
            "log_dir": log_dir.resolve(),
        },
    )

    try:
        if not input_path.exists():
            raise FileNotFoundError(f"Selected input file does not exist: {input_path}")

        df = pd.read_csv(input_path, low_memory=False)
        output = fill_dataset(
            df,
            bottle_url_template=args.bottle_url_template,
            nutrient_source=args.nutrient_source,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output.to_csv(output_path, index=False)
    except Exception:
        LOGGER.exception("Fill failed")
        return 1

    LOGGER.info("Fill completed. Output: %s", output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
