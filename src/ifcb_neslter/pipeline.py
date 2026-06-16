"""End-to-end IFCB processing pipeline."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from .casts import aggregate_cast_data
from .constants import DEFAULT_BOTTLE_URL_TEMPLATE, DEFAULT_TAXONOMY_URL, DEFAULT_STATION_REF_URL
from .metadata import process_meta
from .normalize import normalize
from .taxonomy import import_google_sheet

LOGGER = logging.getLogger("ifcb_neslter")


def validate_required_files(input_dir: Path, required_files: Iterable[str]) -> None:
    """Raise if any required file is missing."""
    missing = [name for name in required_files if not (input_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required file(s) in {input_dir}: {missing}")


def process_dataset(
    input_dir: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    data_type: str,
    meta: pd.DataFrame,
    data_cols: Sequence[str],
    station_reference: str = DEFAULT_STATION_REF_URL,
    max_station_distance_km: float | None = 2.0,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    skip_bottle_merge: bool = False,
) -> Path:
    """Process one raw IFCB dataset and write a cleaned CSV."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    raw_path = input_dir / f"ifcb_{data_type}_raw.csv"
    output_path = output_dir / f"ifcb_{data_type}_clean.csv"

    LOGGER.info("Processing %s data: %s", data_type, raw_path)
    raw = pd.read_csv(raw_path)
    raw_data_cols = [col for col in raw.columns if col in data_cols]

    if "pid" not in meta.columns or "pid" not in raw.columns:
        raise ValueError("Both metadata and raw data must contain a 'pid' column.")

    df = pd.merge(meta, raw, on="pid", how="left")
    df = aggregate_cast_data(
        df,
        raw_data_cols,
        station_reference=station_reference,
        max_station_distance_km=max_station_distance_km,
        bottle_url_template=bottle_url_template,
        skip_bottle_merge=skip_bottle_merge,
    )
    df = normalize(df, raw_data_cols, data_type)

    if "nearest_station" in df.columns:
        df["nearest_station"] = df["nearest_station"].astype("string")
    if "station_distance" in df.columns:
        df["station_distance"] = pd.to_numeric(df["station_distance"], errors="coerce")

    output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    LOGGER.info("Saved cleaned %s data to: %s", data_type, output_path)
    return output_path


def process_all(
    input_dir: str | os.PathLike[str],
    output_dir: str | os.PathLike[str] | None = None,
    sample_type: Sequence[str] | None = None,
    download_taxonomy_if_missing: bool = False,
    taxonomy_url: str = DEFAULT_TAXONOMY_URL,
    station_reference: str = DEFAULT_STATION_REF_URL,
    max_station_distance_km: float | None = 2.0,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    skip_bottle_merge: bool = False,
    data_types: Sequence[str] = ("count", "carbon"),
) -> list[Path]:
    """Run the full IFCB pipeline for selected data types."""
    input_dir = Path(input_dir)
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

    meta = pd.read_csv(input_dir / "ifcb_metadata.csv")
    meta = process_meta(meta, sample_type=sample_type)

    tax = pd.read_csv(taxonomy_path)
    if "Annotations" not in tax.columns:
        raise ValueError("ifcb_taxonomy.csv must contain an 'Annotations' column.")
    data_cols = tax["Annotations"].dropna().astype(str).tolist()

    return [
        process_dataset(
            input_dir=input_dir,
            output_dir=output_dir,
            data_type=data_type,
            meta=meta,
            data_cols=data_cols,
            station_reference=station_reference,
            max_station_distance_km=max_station_distance_km,
            bottle_url_template=bottle_url_template,
            skip_bottle_merge=skip_bottle_merge,
        )
        for data_type in data_types
    ]
