"""Fill missing IFCB cast samples from nearby underway observations."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
import sys
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .constants import DEFAULT_BOTTLE_URL_TEMPLATE, DEFAULT_DATASET
from ifcb.process.logging_utils import log_run_configuration, redact_command_line, setup_logging

LOGGER = logging.getLogger("ifcb.process.neslter")

DEFAULT_NUTRIENT_URL = "https://nes-lter-api.whoi.edu/api/nut/all.csv"
DEFAULT_TARGET_STATIONS = ("MVCO", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11")
DEFAULT_TAXONOMIC_LEVELS = ("Phylum", "Class", "Order", "Family", "Genus", "Species")
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


def add_taxonomy_label(
    taxonomy: pd.DataFrame,
    taxonomic_levels: Sequence[str] = DEFAULT_TAXONOMIC_LEVELS,
    label_col: str = "Label",
) -> pd.DataFrame:
    """Create a filled taxon label from the deepest available taxonomic level.

    For each annotation j, the label is the last non-missing value along:
    Phylum -> Class -> Order -> Family -> Genus -> Species.
    """
    taxonomy = taxonomy.copy()
    for col in taxonomy.columns:
        if pd.api.types.is_object_dtype(taxonomy[col]) or pd.api.types.is_string_dtype(taxonomy[col]):
            taxonomy[col] = taxonomy[col].map(lambda x: x.strip() if isinstance(x, str) else x)
    levels = [col for col in taxonomic_levels if col in taxonomy.columns]
    if not levels:
        raise ValueError("taxonomy must include at least one taxonomic level column.")
    taxonomy[label_col] = taxonomy[levels].ffill(axis=1).iloc[:, -1]
    return taxonomy


def map_taxa_to_label(
    df: pd.DataFrame,
    taxonomy: pd.DataFrame,
    from_level: str = "Annotations",
    to_level: str = "Label",
    aggfunc: str = "sum",
    drop_labels: Iterable[str] = ("Nanoplankton",),
    drop_annotations: Iterable[str] = ("nanoplankton_mix",),
) -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    """Aggregate annotation columns to a taxonomic label while keeping metadata.

    Mathematically, if several annotation columns j map to the same label g,
    the filled table stores X_g = sum_j X_j for each sample row.
    """
    taxonomy = taxonomy.copy()
    for col in taxonomy.columns:
        if pd.api.types.is_object_dtype(taxonomy[col]) or pd.api.types.is_string_dtype(taxonomy[col]):
            taxonomy[col] = taxonomy[col].map(lambda x: x.strip() if isinstance(x, str) else x)
    if to_level not in taxonomy.columns:
        levels = [col for col in DEFAULT_TAXONOMIC_LEVELS if col in taxonomy.columns]
        if not levels:
            raise ValueError("taxonomy must include at least one taxonomic level column.")
        taxonomy[to_level] = taxonomy[levels].ffill(axis=1).iloc[:, -1]
    if from_level not in taxonomy.columns or to_level not in taxonomy.columns:
        raise ValueError(f"taxonomy must include '{from_level}' and '{to_level}'.")

    drop_labels = set(drop_labels)
    drop_annotations = set(drop_annotations)
    mapping = taxonomy.set_index(from_level)[to_level].dropna().to_dict()
    mapped_taxa = [col for col in df.columns if col in mapping and col not in drop_annotations]
    metadata_cols = [col for col in df.columns if col not in mapped_taxa]

    taxa = df.loc[:, mapped_taxa].apply(pd.to_numeric, errors="coerce").copy()
    taxa.columns = [mapping[col] for col in mapped_taxa]
    taxa_agg = taxa.T.groupby(level=0).agg(aggfunc).T
    taxa_agg = taxa_agg.loc[:, [col for col in taxa_agg.columns if col not in drop_labels]]

    out = pd.concat([df.loc[:, metadata_cols].copy(), taxa_agg], axis=1)
    LOGGER.info("Mapped %s annotation columns to %s taxon labels", len(mapped_taxa), len(taxa_agg.columns))
    return out, sorted(taxa_agg.columns.tolist()), taxonomy


def assign_underway_nearest_stations(
    df: pd.DataFrame,
    station_reference: pd.DataFrame | str | os.PathLike[str] | None = None,
    max_station_distance_km: float | None = 2.0,
    sample_type_col: str = "sample_type",
    underway_label: str = "underway",
) -> pd.DataFrame:
    """Assign nearest station to underway rows so they can fill missing casts."""
    from .stations import assign_nearest_stations

    out = df.copy()
    out["sample_time"] = pd.to_datetime(out["sample_time"], errors="coerce", utc=True)
    mask = out[sample_type_col].eq(underway_label)
    if not mask.any():
        return out

    assigned = assign_nearest_stations(
        out.loc[mask],
        station_reference=station_reference,
        lat_col="latitude",
        lon_col="longitude",
        time_col="sample_time",
        max_distance_km=max_station_distance_km,
    )
    out.loc[mask, "nearest_station"] = assigned["nearest_station"].to_numpy()
    out.loc[mask, "station_distance"] = assigned["station_distance"].to_numpy()
    return out


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


def fill_cast_from_udw_from_bottles(
    df: pd.DataFrame,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
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


def merge_nutrients(
    df: pd.DataFrame,
    nutrient_source: str | os.PathLike[str] = DEFAULT_NUTRIENT_URL,
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


def make_filled_dataset(
    input_dir: str | os.PathLike[str],
    output_dir: str | os.PathLike[str] | None = None,
    data_types: Sequence[str] = ("count", "carbon"),
    input_stage: str = "clean",
    output_stage: str = "fill",
    target_stations: Sequence[str] = DEFAULT_TARGET_STATIONS,
    bottle_url_template: str = DEFAULT_BOTTLE_URL_TEMPLATE,
    nutrient_source: str | os.PathLike[str] = DEFAULT_NUTRIENT_URL,
    taxonomy_file: str | os.PathLike[str] = "ifcb_taxonomy.csv",
    input_files: Mapping[str, str | os.PathLike[str]] | None = None,
    output_files: Mapping[str, str | os.PathLike[str]] | None = None,
) -> list[Path]:
    """Create separate *_fill.csv files from clean IFCB files."""
    input_dir = Path(input_dir)
    output_dir = input_dir if output_dir is None else Path(output_dir)
    taxonomy_path = input_dir / taxonomy_file
    if not taxonomy_path.exists():
        raise FileNotFoundError(f"Selected taxonomy file does not exist: {taxonomy_path}")
    taxonomy = pd.read_csv(taxonomy_path, low_memory=False)
    taxonomy = add_taxonomy_label(taxonomy)
    input_files = input_files or {}
    output_files = output_files or {}

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for data_type in data_types:
        input_path = input_dir / input_files.get(data_type, f"ifcb_{data_type}_{input_stage}.csv")
        if not input_path.exists():
            raise FileNotFoundError(f"Selected {data_type} input file does not exist: {input_path}")

        LOGGER.info("Filling missing %s data from %s", data_type, input_path)
        df = pd.read_csv(input_path, low_memory=False)
        input_rows = len(df)
        df, _, taxonomy = map_taxa_to_label(df, taxonomy)
        df = fill_missing_casts_from_underway(
            df,
            target_stations=target_stations,
            mark_filled_col="_fill_created",
        )
        new_rows = df.loc[df["_fill_created"]].drop(columns=["_fill_created"]).copy()
        existing_rows = df.loc[~df["_fill_created"]].drop(columns=["_fill_created"]).copy()
        if not new_rows.empty:
            new_rows = fill_cast_from_udw_from_bottles(new_rows, bottle_url_template=bottle_url_template)
            new_rows = merge_nutrients(new_rows, nutrient_source=nutrient_source)
        df = pd.concat([existing_rows, new_rows], ignore_index=True)

        output_path = output_dir / output_files.get(data_type, f"ifcb_{data_type}_{output_stage}.csv")
        df.sort_values("sample_time").reset_index(drop=True).to_csv(output_path, index=False)
        outputs.append(output_path)
        LOGGER.info("Saved fill %s data to: %s (%s -> %s rows)", data_type, output_path, input_rows, len(df))

    return outputs


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create separate IFCB *_fill.csv files from existing *_clean.csv files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_data_path", nargs="?", default=None, help="Directory containing IFCB *_clean.csv files.")
    parser.add_argument("--dataset", default=DEFAULT_DATASET, help="Dataset folder under current-directory data/.")
    parser.add_argument("-o", "--output-dir", default=None, help="Directory for fill output CSVs.")
    parser.add_argument("--data-type", choices=["count", "carbon"], nargs="+", default=["count", "carbon"])
    parser.add_argument("--input-stage", default="clean")
    parser.add_argument("--output-stage", default="fill")
    parser.add_argument("--taxonomy-file", default="ifcb_taxonomy.csv", help="Taxonomy input filename.")
    parser.add_argument("--count-file", default=None, help="Count input filename; defaults from --input-stage.")
    parser.add_argument("--carbon-file", default=None, help="Carbon input filename; defaults from --input-stage.")
    parser.add_argument("--count-output-file", default=None, help="Count output filename; defaults from --output-stage.")
    parser.add_argument("--carbon-output-file", default=None, help="Carbon output filename; defaults from --output-stage.")
    parser.add_argument("--bottle-url-template", default=DEFAULT_BOTTLE_URL_TEMPLATE)
    parser.add_argument("--nutrient-source", default=DEFAULT_NUTRIENT_URL)
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    parser.add_argument("--log-dir", default=None, help="Directory for timestamped workflow log files.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.input_data_path is None:
        input_dir = Path.cwd().resolve() / "data" / args.dataset
    else:
        input_dir = Path(args.input_data_path)
    output_dir = Path(args.output_dir) if args.output_dir is not None else input_dir
    log_dir = Path(args.log_dir) if args.log_dir is not None else output_dir / "logs"
    setup_logging(log_dir=log_dir, name="ifcb_fill_missing", level=getattr(logging, args.log_level))
    LOGGER.info("Starting IFCB missing-cast fill for dataset %s", args.dataset)
    log_run_configuration(
        LOGGER,
        {
            "command": redact_command_line(sys.argv),
            "dataset": args.dataset,
            "input_dir": input_dir.resolve(),
            "output_dir": output_dir.resolve(),
            "data_type": args.data_type,
            "input_stage": args.input_stage,
            "output_stage": args.output_stage,
            "taxonomy_file": args.taxonomy_file,
            "count_file": args.count_file,
            "carbon_file": args.carbon_file,
            "count_output_file": args.count_output_file,
            "carbon_output_file": args.carbon_output_file,
            "bottle_url_template": args.bottle_url_template,
            "nutrient_source": args.nutrient_source,
            "log_level": args.log_level,
            "log_dir": log_dir.resolve(),
        },
    )

    try:
        outputs = make_filled_dataset(
            input_dir=input_dir,
            output_dir=args.output_dir,
            data_types=args.data_type,
            input_stage=args.input_stage,
            output_stage=args.output_stage,
            bottle_url_template=args.bottle_url_template,
            nutrient_source=args.nutrient_source,
            taxonomy_file=args.taxonomy_file,
            input_files={
                key: value
                for key, value in {"count": args.count_file, "carbon": args.carbon_file}.items()
                if value is not None
            },
            output_files={
                key: value
                for key, value in {
                    "count": args.count_output_file,
                    "carbon": args.carbon_output_file,
                }.items()
                if value is not None
            },
        )
    except Exception:
        LOGGER.exception("Fill failed")
        return 1

    LOGGER.info("Fill completed. Outputs: %s", ", ".join(str(path) for path in outputs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
