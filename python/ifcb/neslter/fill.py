"""Fill missing IFCB cast samples from nearby underway observations."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .constants import DEFAULT_BOTTLE_URL_TEMPLATE, DEFAULT_DATASET

LOGGER = logging.getLogger("ifcb.neslter")

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


def _read_csv(path_or_url: str | os.PathLike[str]) -> pd.DataFrame:
    """Read local or remote CSV with Pandas' mixed-type warning avoided."""
    return pd.read_csv(path_or_url, low_memory=False)


def _repo_root() -> Path:
    """Return the repository root for an editable source checkout."""
    return Path(__file__).resolve().parents[3]


def _default_data_dir(dataset: str = DEFAULT_DATASET) -> Path:
    """Return the repo-local data directory without importing the process module."""
    return _repo_root() / "data" / dataset


def _strip_strings(df: pd.DataFrame) -> pd.DataFrame:
    """Trim whitespace in string cells while preserving missing values."""
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_object_dtype(out[col]) or pd.api.types.is_string_dtype(out[col]):
            out[col] = out[col].map(lambda x: x.strip() if isinstance(x, str) else x)
    return out


def add_taxonomy_label(
    taxonomy: pd.DataFrame,
    taxonomic_levels: Sequence[str] = DEFAULT_TAXONOMIC_LEVELS,
    label_col: str = "Label",
) -> pd.DataFrame:
    """Create a filled taxon label from the deepest available taxonomic level.

    For each annotation j, the label is the last non-missing value along:
    Phylum -> Class -> Order -> Family -> Genus -> Species.
    """
    taxonomy = _strip_strings(taxonomy)
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
    taxonomy = add_taxonomy_label(taxonomy) if to_level not in taxonomy.columns else _strip_strings(taxonomy)
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
) -> pd.DataFrame:
    """Create synthetic surface casts from underway rows when station casts are missing.

    Ecological goal: preserve a balanced local-community set i across surveys by
    filling absent surface observations with the nearest underway sample from
    the same cruise and station.
    """
    df = df.copy()
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
        filled.append(best)

    if filled:
        df = pd.concat([df, pd.DataFrame(filled)], ignore_index=True)
    return df.sort_values([cruise_col, station_col, cast_col, time_col, depth_col]).reset_index(drop=True)


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
    out["_cruise_key"] = out[cruise_col].astype(str).str.strip()
    out["_cast_key"] = pd.to_numeric(out[cast_col], errors="coerce").astype("Int64")
    out["_sample_type_key"] = out[sample_type_col].astype(str).str.strip().str.lower()
    target_mask = out["_sample_type_key"].eq(target_sample_type.lower()) & out["_cast_key"].notna()

    for cruise in out.loc[target_mask, "_cruise_key"].dropna().unique():
        url = bottle_url_template.format(cruise=str(cruise).lower())
        try:
            bottle = _read_csv(url)
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

    return out.drop(columns=["_cruise_key", "_cast_key", "_sample_type_key"], errors="ignore")


def merge_nutrients(
    df: pd.DataFrame,
    nutrient_source: str | os.PathLike[str] = DEFAULT_NUTRIENT_URL,
    nutrient_cols: Sequence[str] = DEFAULT_NUTRIENT_COLS,
    cruise_col: str = "cruise",
    cast_col: str = "cast",
    niskin_col: str = "niskin",
) -> pd.DataFrame:
    """Attach nutrient values by cruise, cast, and best matching niskin bottle."""
    nut = _read_csv(nutrient_source)
    available_cols = [col for col in nutrient_cols if col in nut.columns]
    if not available_cols:
        LOGGER.warning("No expected nutrient columns found in %s.", nutrient_source)
        return df.copy()

    nut_avg = nut.groupby([cruise_col, niskin_col, cast_col], as_index=False)[available_cols].mean()
    nut_avg[cast_col] = pd.to_numeric(nut_avg[cast_col], errors="coerce")
    nut_avg[niskin_col] = pd.to_numeric(nut_avg[niskin_col], errors="coerce")

    flag_cols = [col for col in available_cols if col.startswith("flag_")]
    nut_avg["_flag_score"] = nut_avg[flag_cols].sum(axis=1, skipna=True) if flag_cols else 0
    if flag_cols:
        nut_avg.loc[nut_avg[flag_cols].notna().sum(axis=1).eq(0), "_flag_score"] = np.inf

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
    merged["_matched"] = merged[list(source_cols.values())].notna().any(axis=1).astype(int)
    merged = merged.sort_values(["_row_id", "_matched", "_flag_score"], ascending=[True, False, True])
    best = merged.groupby("_row_id", as_index=False).first()

    for col in available_cols:
        source_col = source_cols[col]
        if source_col in best.columns:
            out[col] = best[source_col].to_numpy()

    matched = best["_matched"].eq(1).to_numpy()
    if f"{niskin_col}_nut" in best.columns:
        out.loc[matched, niskin_col] = best.loc[matched, f"{niskin_col}_nut"].astype("Int64").astype(str).to_numpy()
    elif niskin_col in best.columns:
        out.loc[matched, niskin_col] = best.loc[matched, niskin_col].astype("Int64").astype(str).to_numpy()

    return out.drop(columns=["_row_id", "_cast_numeric", "_niskin_candidates"], errors="ignore")


def make_filled_dataset(
    input_dir: str | os.PathLike[str],
    output_dir: str | os.PathLike[str] | None = None,
    data_types: Sequence[str] = ("carbon",),
    input_stage: str = "station",
    output_stage: str = "filled",
    target_stations: Sequence[str] = DEFAULT_TARGET_STATIONS,
) -> list[Path]:
    """Create separate *_filled.csv files from station-assigned IFCB files."""
    input_dir = Path(input_dir)
    output_dir = input_dir if output_dir is None else Path(output_dir)
    taxonomy = _read_csv(input_dir / "ifcb_taxonomy.csv")
    taxonomy = add_taxonomy_label(taxonomy)

    output_dir.mkdir(parents=True, exist_ok=True)
    taxonomy_path = output_dir / "ifcb_taxonomy_filled.csv"
    taxonomy.to_csv(taxonomy_path, index=False)

    outputs = [taxonomy_path]
    for data_type in data_types:
        input_path = input_dir / f"ifcb_{data_type}_{input_stage}.csv"
        if not input_path.exists():
            raise FileNotFoundError(f"Missing station-assigned input file: {input_path}")

        LOGGER.info("Filling missing %s data from %s", data_type, input_path)
        df = _read_csv(input_path)
        df, _, taxonomy = map_taxa_to_label(df, taxonomy)
        df = fill_missing_casts_from_underway(df, target_stations=target_stations)

        output_path = output_dir / f"ifcb_{data_type}_{output_stage}.csv"
        df.sort_values("sample_time").reset_index(drop=True).to_csv(output_path, index=False)
        outputs.append(output_path)
        LOGGER.info("Saved filled %s data to: %s", data_type, output_path)

    return outputs


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create separate IFCB *_filled.csv files from existing *_clean.csv files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_data_path", nargs="?", default=None, help="Directory containing station-assigned IFCB CSV files.")
    parser.add_argument("--dataset", default=DEFAULT_DATASET, help="Dataset folder under repo-local data/.")
    parser.add_argument("-o", "--output-dir", default=None, help="Directory for filled output CSVs.")
    parser.add_argument("--data-type", choices=["count", "carbon"], nargs="+", default=["carbon"])
    parser.add_argument("--input-stage", default="station")
    parser.add_argument("--output-stage", default="filled")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    input_dir = _default_data_dir(args.dataset) if args.input_data_path is None else Path(args.input_data_path)

    try:
        outputs = make_filled_dataset(
            input_dir=input_dir,
            output_dir=args.output_dir,
            data_types=args.data_type,
            input_stage=args.input_stage,
            output_stage=args.output_stage,
        )
    except Exception as exc:
        LOGGER.error("Fill failed: %s", exc)
        return 1

    LOGGER.info("Fill completed. Outputs: %s", ", ".join(str(path) for path in outputs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
