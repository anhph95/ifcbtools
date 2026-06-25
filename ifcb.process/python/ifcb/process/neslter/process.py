"""Data processing for MATLAB-exported IFCB CSV files."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Sequence

import pandas as pd

from .casts import aggregate_cast_data
from .constants import (
    DEFAULT_BOTTLE_URL_TEMPLATE,
    DEFAULT_STATION_REF_URL,
    DEFAULT_TAXONOMY_URL,
)
from .fill import DEFAULT_NUTRIENT_URL, merge_nutrients
from .metadata import process_meta
from .normalize import normalize
from .pipeline import add_nearest_station, merge_bottle
from .taxonomy import import_google_sheet

LOGGER = logging.getLogger("ifcb.process.neslter")


def normalization_scaling_factor(input_file: str | os.PathLike[str]) -> float:
    """Return the unit conversion scale implied by the IFCB product filename."""
    stem = Path(input_file).stem.lower()
    if "count" in stem:
        return 1000.0
    if "carbon" in stem:
        return 0.001
    raise ValueError(
        "Cannot infer IFCB normalization scale from input filename. "
        "Use a count or carbon filename such as ifcb_count.csv or ifcb_carbon.csv."
    )


def default_output_path(
    input_file: str | os.PathLike[str],
    *,
    clean: bool = False,
    add_station: bool = False,
    merge_bottle_data: bool = False,
    merge_nutrient: bool = False,
) -> Path:
    """Append selected operation suffixes to an input filename."""
    input_path = Path(input_file)
    suffixes = []
    if clean:
        suffixes.append("clean")
    if add_station:
        suffixes.append("station")
    if merge_bottle_data:
        suffixes.append("bottle")
    if merge_nutrient:
        suffixes.append("nutrient")
    operation_suffix = "".join(f"_{suffix}" for suffix in suffixes)
    return input_path.with_name(f"{input_path.stem}{operation_suffix}{input_path.suffix}")


def process(
    input_file: str | os.PathLike[str],
    output_file: str | os.PathLike[str] | None = None,
    sample_type: Sequence[str] | None = None,
    download_taxonomy_if_missing: bool = True,
    taxonomy_url: str = DEFAULT_TAXONOMY_URL,
    station_reference: str | os.PathLike[str] | None = DEFAULT_STATION_REF_URL,
    max_station_distance_km: float | None = 2.0,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    nutrient_source: str | os.PathLike[str] = DEFAULT_NUTRIENT_URL,
    metadata_file: str | os.PathLike[str] | None = None,
    taxonomy_file: str | os.PathLike[str] | None = None,
    clean: bool = False,
    add_station: bool = False,
    merge_bottle_data: bool = False,
    merge_nutrient: bool = False,
) -> Path:
    """Run selected cleaning and enrichment operations on one IFCB CSV file."""
    input_path = Path(input_file)
    output_path = (
        Path(output_file)
        if output_file is not None
        else default_output_path(
            input_path,
            clean=clean,
            add_station=add_station,
            merge_bottle_data=merge_bottle_data,
            merge_nutrient=merge_nutrient,
        )
    )
    metadata_path = Path(metadata_file) if metadata_file is not None else input_path.parent / "ifcb_metadata.csv"
    taxonomy_path = Path(taxonomy_file) if taxonomy_file is not None else input_path.parent / "ifcb_taxonomy.csv"

    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")

    if not any((clean, add_station, merge_bottle_data, merge_nutrient)):
        raise ValueError("At least one processing operation must be selected.")

    if clean:
        if not taxonomy_path.exists() and download_taxonomy_if_missing:
            import_google_sheet(taxonomy_url, save_path=taxonomy_path)

        selected_inputs = [metadata_path, taxonomy_path]
        missing = [path for path in selected_inputs if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Selected input file(s) do not exist: {missing}")

        LOGGER.info("Reading metadata: %s", metadata_path)
        meta = pd.read_csv(metadata_path, low_memory=False)
        meta = process_meta(meta, sample_type=sample_type)
        LOGGER.info("Prepared metadata rows: %s", len(meta))

        tax = pd.read_csv(taxonomy_path)
        if "Annotations" not in tax.columns:
            raise ValueError(f"{taxonomy_path} must contain an 'Annotations' column.")
        data_cols = tax["Annotations"].dropna().astype(str).tolist()
        LOGGER.info("Loaded taxonomy annotations: %s", len(data_cols))

        LOGGER.info("Cleaning data: %s", input_path)
        raw = pd.read_csv(input_path)
        LOGGER.info("Read %s raw rows and %s metadata rows", len(raw), len(meta))
        raw_data_cols = [col for col in raw.columns if col in data_cols]
        LOGGER.info("Detected %s taxon columns", len(raw_data_cols))

        if "pid" not in meta.columns or "pid" not in raw.columns:
            raise ValueError("Both metadata and raw data must contain a 'pid' column.")

        df = pd.merge(meta, raw, on="pid", how="left")
        LOGGER.info("Merged metadata and data: %s rows", len(df))
        df = aggregate_cast_data(df, raw_data_cols)
        scaling_factor = normalization_scaling_factor(input_path)
        LOGGER.info("Normalizing taxon columns with scaling factor: %s", scaling_factor)
        df = normalize(df, raw_data_cols, scaling_factor=scaling_factor)
    else:
        LOGGER.info("Reading input: %s", input_path)
        df = pd.read_csv(input_path, low_memory=False)

    if add_station:
        df = add_nearest_station(
            df,
            station_reference=station_reference,
            max_station_distance_km=max_station_distance_km,
        )
    if merge_bottle_data:
        df = merge_bottle(df, bottle_url_template=bottle_url_template)
    if merge_nutrient:
        df = merge_nutrients(df, nutrient_source=nutrient_source)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    LOGGER.info("Saved processed data to: %s", output_path)
    return output_path
