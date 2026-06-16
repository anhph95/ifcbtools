"""Explicit post-clean IFCB processing stages."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Sequence

import pandas as pd

from .casts import normalize_integer_key
from .constants import DEFAULT_BOTTLE_URL_TEMPLATE, DEFAULT_DATASET
from .fill import DEFAULT_NUTRIENT_URL, fill_cast_from_udw_from_bottles, merge_nutrients

LOGGER = logging.getLogger("ifcb.neslter")


def repo_root() -> Path:
    """Return the repository root for an editable source checkout."""
    return Path(__file__).resolve().parents[3]


def default_data_dir(dataset: str = DEFAULT_DATASET) -> Path:
    """Return the repo-local data directory for a dataset."""
    return repo_root() / "data" / dataset


def data_path(data_dir: str | os.PathLike[str], data_type: str, stage: str) -> Path:
    """Build a staged IFCB file path such as ifcb_carbon_station.csv."""
    return Path(data_dir) / f"ifcb_{data_type}_{stage}.csv"


def add_nearest_station(
    df: pd.DataFrame,
    station_reference: pd.DataFrame | str | os.PathLike[str] | None = None,
    max_station_distance_km: float | None = 2.0,
) -> pd.DataFrame:
    """Assign nearest_station and station_distance to every row using StationLocator."""
    from .stations import StationLocator

    locator = StationLocator(station_reference=station_reference, max_distance_km=max_station_distance_km)
    out = df.copy()
    out["sample_time"] = pd.to_datetime(out["sample_time"], errors="coerce", utc=True)

    names, distances = locator.nearest_stations(
        out,
        lat_col="latitude",
        lon_col="longitude",
        time_col="sample_time",
        max_distance_km=max_station_distance_km,
    )
    out["nearest_station"] = names
    out["station_distance"] = distances
    return out


def add_nearest_station_files(
    data_dir: str | os.PathLike[str],
    output_dir: str | os.PathLike[str] | None = None,
    data_types: Sequence[str] = ("count", "carbon"),
    input_stage: str = "clean",
    output_stage: str = "station",
    station_reference: pd.DataFrame | str | os.PathLike[str] | None = None,
    max_station_distance_km: float | None = 2.0,
) -> list[Path]:
    """Add nearest-station columns to staged IFCB files."""
    data_dir = Path(data_dir)
    output_dir = data_dir if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    for data_type in data_types:
        input_path = data_path(data_dir, data_type, input_stage)
        output_path = data_path(output_dir, data_type, output_stage)
        df = pd.read_csv(input_path, low_memory=False)
        add_nearest_station(
            df,
            station_reference=station_reference,
            max_station_distance_km=max_station_distance_km,
        ).to_csv(output_path, index=False)
        outputs.append(output_path)
        LOGGER.info("Saved station-assigned %s data to: %s", data_type, output_path)

    return outputs


def merge_bottle(
    df: pd.DataFrame,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    cruise_col: str = "cruise",
    cast_col: str = "cast",
    niskin_col: str = "niskin",
) -> pd.DataFrame:
    """Merge CTD bottle data by cruise, cast, and niskin where possible."""
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
        for key in [cast_col, niskin_col]:
            if key in bottle.columns:
                bottle[key] = normalize_integer_key(bottle[key])
            if key in sub.columns:
                sub[key] = normalize_integer_key(sub[key])

        if cast_col not in bottle.columns or niskin_col not in bottle.columns:
            LOGGER.warning("Bottle data for %s lacks cast/niskin keys; keeping rows unchanged.", cruise)
            out_frames.append(sub)
            continue

        bottle = bottle.drop_duplicates(subset=[cast_col, niskin_col])
        merged = pd.merge(
            sub,
            bottle,
            on=[cast_col, niskin_col],
            how="left",
            suffixes=("", "_bottle"),
            validate="many_to_one",
        )
        out_frames.append(merged)

    out = pd.concat(out_frames, ignore_index=True).sort_values("sample_time").reset_index(drop=True)
    return fill_cast_from_udw_from_bottles(out, bottle_url_template=bottle_url_template)


def merge_bottle_files(
    data_dir: str | os.PathLike[str],
    output_dir: str | os.PathLike[str] | None = None,
    data_types: Sequence[str] = ("count", "carbon"),
    input_stage: str = "filled",
    output_stage: str = "bottle",
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
) -> list[Path]:
    """Merge CTD bottle data into staged IFCB files."""
    data_dir = Path(data_dir)
    output_dir = data_dir if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    for data_type in data_types:
        input_path = data_path(data_dir, data_type, input_stage)
        output_path = data_path(output_dir, data_type, output_stage)
        df = pd.read_csv(input_path, low_memory=False)
        merge_bottle(df, bottle_url_template=bottle_url_template).to_csv(output_path, index=False)
        outputs.append(output_path)
        LOGGER.info("Saved bottle-merged %s data to: %s", data_type, output_path)

    return outputs


def merge_nutrient_files(
    data_dir: str | os.PathLike[str],
    output_dir: str | os.PathLike[str] | None = None,
    data_types: Sequence[str] = ("carbon",),
    input_stage: str = "bottle",
    output_stage: str = "nutrient",
    nutrient_source: str | os.PathLike[str] = DEFAULT_NUTRIENT_URL,
) -> list[Path]:
    """Merge nutrient data into staged IFCB files."""
    data_dir = Path(data_dir)
    output_dir = data_dir if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    for data_type in data_types:
        input_path = data_path(data_dir, data_type, input_stage)
        output_path = data_path(output_dir, data_type, output_stage)
        df = pd.read_csv(input_path, low_memory=False)
        merge_nutrients(df, nutrient_source=nutrient_source).to_csv(output_path, index=False)
        outputs.append(output_path)
        LOGGER.info("Saved nutrient-merged %s data to: %s", data_type, output_path)

    return outputs


def _common_stage_parser(description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_data_path", nargs="?", default=None, help="Directory containing staged IFCB CSV files.")
    parser.add_argument("--dataset", default=DEFAULT_DATASET, help="Dataset folder under repo-local data/.")
    parser.add_argument("-o", "--output-dir", default=None, help="Directory for output CSVs.")
    parser.add_argument("--data-type", choices=["count", "carbon"], nargs="+", default=["count", "carbon"])
    parser.add_argument("--input-stage", default=None)
    parser.add_argument("--output-stage", default=None)
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    return parser


def add_stations_main(argv: Sequence[str] | None = None) -> int:
    parser = _common_stage_parser("Add nearest_station and station_distance to all IFCB rows.")
    parser.set_defaults(input_stage="clean", output_stage="station")
    parser.add_argument("--station-reference", default=None, help="Station reference CSV path.")
    parser.add_argument("--max-station-distance-km", type=float, default=2.0)
    parser.add_argument("--no-station-distance-limit", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    data_dir = default_data_dir(args.dataset) if args.input_data_path is None else Path(args.input_data_path)
    max_distance = None if args.no_station_distance_limit else args.max_station_distance_km
    input_stage = args.input_stage or "clean"
    output_stage = args.output_stage or "station"
    try:
        outputs = add_nearest_station_files(
            data_dir=data_dir,
            output_dir=args.output_dir,
            data_types=args.data_type,
            input_stage=input_stage,
            output_stage=output_stage,
            station_reference=args.station_reference,
            max_station_distance_km=max_distance,
        )
    except Exception as exc:
        LOGGER.error("Add-stations failed: %s", exc)
        return 1
    LOGGER.info("Add-stations completed. Outputs: %s", ", ".join(str(path) for path in outputs))
    return 0


def merge_bottle_main(argv: Sequence[str] | None = None) -> int:
    parser = _common_stage_parser("Merge CTD bottle data into staged IFCB files.")
    parser.set_defaults(input_stage="filled", output_stage="bottle")
    parser.add_argument("--bottle-url-template", default=DEFAULT_BOTTLE_URL_TEMPLATE)
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    data_dir = default_data_dir(args.dataset) if args.input_data_path is None else Path(args.input_data_path)
    input_stage = args.input_stage or "filled"
    output_stage = args.output_stage or "bottle"
    try:
        outputs = merge_bottle_files(
            data_dir=data_dir,
            output_dir=args.output_dir,
            data_types=args.data_type,
            input_stage=input_stage,
            output_stage=output_stage,
            bottle_url_template=args.bottle_url_template,
        )
    except Exception as exc:
        LOGGER.error("Merge-bottle failed: %s", exc)
        return 1
    LOGGER.info("Merge-bottle completed. Outputs: %s", ", ".join(str(path) for path in outputs))
    return 0


def merge_nutrient_main(argv: Sequence[str] | None = None) -> int:
    parser = _common_stage_parser("Merge nutrient data into staged IFCB files.")
    parser.set_defaults(input_stage="bottle", output_stage="nutrient")
    parser.add_argument("--nutrient-source", default=DEFAULT_NUTRIENT_URL)
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    data_dir = default_data_dir(args.dataset) if args.input_data_path is None else Path(args.input_data_path)
    input_stage = args.input_stage or "bottle"
    output_stage = args.output_stage or "nutrient"
    try:
        outputs = merge_nutrient_files(
            data_dir=data_dir,
            output_dir=args.output_dir,
            data_types=args.data_type,
            input_stage=input_stage,
            output_stage=output_stage,
            nutrient_source=args.nutrient_source,
        )
    except Exception as exc:
        LOGGER.error("Merge-nutrient failed: %s", exc)
        return 1
    LOGGER.info("Merge-nutrient completed. Outputs: %s", ", ".join(str(path) for path in outputs))
    return 0
