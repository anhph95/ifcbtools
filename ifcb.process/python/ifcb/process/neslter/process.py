"""Data processing for MATLAB-exported IFCB CSV files."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from .casts import aggregate_cast_data
from .constants import (
    DEFAULT_BOTTLE_URL_TEMPLATE,
    DEFAULT_DATASET,
    DEFAULT_STATION_REF_URL,
    DEFAULT_TAXONOMY_URL,
)
from .fill import DEFAULT_NUTRIENT_URL, merge_nutrients
from .metadata import process_meta
from .normalize import normalize
from .pipeline import add_nearest_station, merge_bottle
from .taxonomy import import_google_sheet

LOGGER = logging.getLogger("ifcb.process.neslter")


def repo_root() -> Path:
    """Return the repository root for an editable source checkout."""
    return Path(__file__).resolve().parents[5]


def matlab_export_data_dir(dataset: str = DEFAULT_DATASET) -> Path:
    """Return the default data directory written by ifcb.process/matlab/export_ifcb_mat.m."""
    return repo_root() / "data" / dataset


def validate_required_files(input_dir: Path, required_files: Iterable[str]) -> None:
    """Raise if any required file is missing."""
    missing = [name for name in required_files if not (input_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required file(s) in {input_dir}: {missing}")


def process_data_type(
    input_dir: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    data_type: str,
    meta: pd.DataFrame,
    data_cols: Sequence[str],
    station_reference: str | os.PathLike[str] | None = DEFAULT_STATION_REF_URL,
    max_station_distance_km: float | None = 2.0,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    nutrient_source: str | os.PathLike[str] = DEFAULT_NUTRIENT_URL,
) -> Path:
    """Process one raw IFCB data type through the clean pipeline and write CSV."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_path = output_dir / f"ifcb_{data_type}_clean.csv"

    raw_path = input_dir / f"ifcb_{data_type}_raw.csv"
    LOGGER.info("Processing %s data: %s", data_type, raw_path)
    raw = pd.read_csv(raw_path)
    LOGGER.info("Read %s raw rows and %s metadata rows", len(raw), len(meta))
    raw_data_cols = [col for col in raw.columns if col in data_cols]
    LOGGER.info("Detected %s %s taxon columns", len(raw_data_cols), data_type)

    if "pid" not in meta.columns or "pid" not in raw.columns:
        raise ValueError("Both metadata and raw data must contain a 'pid' column.")

    df = pd.merge(meta, raw, on="pid", how="left")
    LOGGER.info("Merged metadata and %s data: %s rows", data_type, len(df))
    df = aggregate_cast_data(df, raw_data_cols)
    df = normalize(df, raw_data_cols, data_type)
    df = add_nearest_station(
        df,
        station_reference=station_reference,
        max_station_distance_km=max_station_distance_km,
    )
    df = merge_bottle(df, bottle_url_template=bottle_url_template)
    df = merge_nutrients(df, nutrient_source=nutrient_source)

    output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    LOGGER.info("Saved cleaned %s data to: %s", data_type, output_path)
    return output_path


def process(
    input_dir: str | os.PathLike[str] | None = None,
    output_dir: str | os.PathLike[str] | None = None,
    dataset: str = DEFAULT_DATASET,
    sample_type: Sequence[str] | None = None,
    download_taxonomy_if_missing: bool = True,
    taxonomy_url: str = DEFAULT_TAXONOMY_URL,
    station_reference: str | os.PathLike[str] | None = DEFAULT_STATION_REF_URL,
    max_station_distance_km: float | None = 2.0,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    nutrient_source: str | os.PathLike[str] = DEFAULT_NUTRIENT_URL,
    data_types: Sequence[str] = ("count", "carbon"),
) -> list[Path]:
    """Process MATLAB-exported IFCB CSV files for selected data types."""
    input_dir = matlab_export_data_dir(dataset) if input_dir is None else Path(input_dir)
    output_dir = Path(output_dir) if output_dir is not None else input_dir
    taxonomy_path = input_dir / "ifcb_taxonomy.csv"

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    if not taxonomy_path.exists() and download_taxonomy_if_missing:
        import_google_sheet(taxonomy_url, save_path=taxonomy_path)

    required_files = ["ifcb_metadata.csv", "ifcb_taxonomy.csv"] + [
        f"ifcb_{data_type}_raw.csv" for data_type in data_types
    ]
    validate_required_files(input_dir, required_files)

    LOGGER.info("Reading metadata: %s", input_dir / "ifcb_metadata.csv")
    meta = pd.read_csv(input_dir / "ifcb_metadata.csv", low_memory=False)
    meta = process_meta(meta, sample_type=sample_type)
    LOGGER.info("Prepared metadata rows: %s", len(meta))

    tax = pd.read_csv(taxonomy_path)
    if "Annotations" not in tax.columns:
        raise ValueError("ifcb_taxonomy.csv must contain an 'Annotations' column.")
    data_cols = tax["Annotations"].dropna().astype(str).tolist()
    LOGGER.info("Loaded taxonomy annotations: %s", len(data_cols))

    outputs: list[Path] = []
    station_lookup: pd.DataFrame | None = None
    for data_type in data_types:
        output_path = output_dir / f"ifcb_{data_type}_clean.csv"
        raw_path = input_dir / f"ifcb_{data_type}_raw.csv"
        LOGGER.info("Processing %s data: %s", data_type, raw_path)
        raw = pd.read_csv(raw_path)
        LOGGER.info("Read %s raw rows and %s metadata rows", len(raw), len(meta))
        raw_data_cols = [col for col in raw.columns if col in data_cols]
        LOGGER.info("Detected %s %s taxon columns", len(raw_data_cols), data_type)

        if "pid" not in meta.columns or "pid" not in raw.columns:
            raise ValueError("Both metadata and raw data must contain a 'pid' column.")

        df = pd.merge(meta, raw, on="pid", how="left")
        LOGGER.info("Merged metadata and %s data: %s rows", data_type, len(df))
        df = aggregate_cast_data(df, raw_data_cols)
        df = normalize(df, raw_data_cols, data_type)

        if station_lookup is not None and "pid" in df.columns and set(df["pid"]).issubset(set(station_lookup["pid"])):
            df = df.drop(columns=["nearest_station", "station_distance"], errors="ignore").merge(
                station_lookup,
                on="pid",
                how="left",
                validate="many_to_one",
            )
            LOGGER.info("Reused station assignments for %s by pid: %s rows", data_type, df["nearest_station"].notna().sum())
        else:
            df = add_nearest_station(
                df,
                station_reference=station_reference,
                max_station_distance_km=max_station_distance_km,
            )
            if "pid" in df.columns:
                station_lookup = df[["pid", "nearest_station", "station_distance"]].drop_duplicates("pid")
                LOGGER.info("Built station lookup by pid: %s rows", len(station_lookup))

        df = merge_bottle(df, bottle_url_template=bottle_url_template)
        df = merge_nutrients(df, nutrient_source=nutrient_source)
        output_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)
        LOGGER.info("Saved cleaned %s data to: %s", data_type, output_path)
        outputs.append(output_path)

    return outputs
