"""Data processing for MATLAB-exported IFCB CSV files."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from .casts import aggregate_cast_data
from .constants import (
    DEFAULT_DATASET,
    DEFAULT_TAXONOMY_URL,
)
from .metadata import process_meta
from .normalize import normalize
from .taxonomy import import_google_sheet

LOGGER = logging.getLogger("ifcb.neslter")


def repo_root() -> Path:
    """Return the repository root for an editable source checkout."""
    return Path(__file__).resolve().parents[3]


def matlab_export_data_dir(dataset: str = DEFAULT_DATASET) -> Path:
    """Return the default data directory written by matlab/export_ifcb_mat.m."""
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
) -> Path:
    """Process one raw IFCB data type and write a cleaned CSV."""
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
    df = aggregate_cast_data(df, raw_data_cols)
    df = normalize(df, raw_data_cols, data_type)

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

    meta = pd.read_csv(input_dir / "ifcb_metadata.csv")
    meta = process_meta(meta, sample_type=sample_type)

    tax = pd.read_csv(taxonomy_path)
    if "Annotations" not in tax.columns:
        raise ValueError("ifcb_taxonomy.csv must contain an 'Annotations' column.")
    data_cols = tax["Annotations"].dropna().astype(str).tolist()

    return [
        process_data_type(
            input_dir=input_dir,
            output_dir=output_dir,
            data_type=data_type,
            meta=meta,
            data_cols=data_cols,
        )
        for data_type in data_types
    ]
