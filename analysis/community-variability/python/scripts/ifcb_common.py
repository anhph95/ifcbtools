"""Shared IFCB community-variability workflow steps."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[4]
COMMUNITY_VARIABILITY_SRC = REPO_ROOT / "community.variability" / "python"
IFCB_PROCESS_SRC = REPO_ROOT / "ifcb.process" / "python"
if str(COMMUNITY_VARIABILITY_SRC) not in sys.path:
    sys.path.insert(0, str(COMMUNITY_VARIABILITY_SRC))
if str(IFCB_PROCESS_SRC) not in sys.path:
    sys.path.insert(0, str(IFCB_PROCESS_SRC))

from community_variability import make_community_array
from ifcb.process.logging_utils import log_run_configuration, redact_command_line, setup_logging


SEASONS = ["JFM", "AMJ", "JAS", "OND"]
STATION_LIST = ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11"]
MAIN_CRUISE = [
    "EN608",
    "AR28B",
    "EN617",
    "AR32",
    "EN627",
    "AR34B",
    "EN644",
    "AR39B",
    "EN649",
    "AR44",
    "EN655",
    "EN657",
    "EN661",
    "AR52B",
    "EN668",
    "AR61B",
    "AT46",
    "AR66B",
    "EN687",
    "AR70B",
    "EN695",
    "HRS2303",
    "EN706",
    "AR77",
    "EN712",
    "EN715",
    "EN720",
    "AE2426",
    "EN727",
    "AR88",
    "AR92",
    "AR95",
    "AR99",
]


def repo_root() -> Path:
    """Return the repository root from analysis/community-variability/python/scripts."""
    return REPO_ROOT


def default_data_dir() -> Path:
    """Use the repo-local transect data produced by the Python processing step."""
    return repo_root() / "data" / "NESLTER_transect"


def default_results_dir() -> Path:
    """Use the same results directory as the R analysis project."""
    return repo_root() / "analysis" / "community-variability" / "results"


def add_logging_arguments(parser: argparse.ArgumentParser) -> None:
    """Add consistent logging controls to an analysis workflow parser."""
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    parser.add_argument("--log-dir", default=None, help="Directory for timestamped workflow log files.")


def configure_workflow_logging(args: argparse.Namespace, name: str) -> logging.Logger:
    """Write workflow logs beside analysis results unless explicitly overridden."""
    results_dir = Path(args.results_dir)
    log_dir = Path(args.log_dir) if args.log_dir is not None else results_dir / "logs"
    logger = setup_logging(log_dir=log_dir, name=name, level=getattr(logging, args.log_level))
    settings = vars(args).copy()
    settings["command"] = redact_command_line(sys.argv)
    settings["results_dir"] = results_dir.resolve()
    if "data_dir" in settings:
        settings["data_dir"] = Path(settings["data_dir"]).resolve()
    settings["log_dir"] = log_dir.resolve()
    log_run_configuration(logger, settings)
    return logger


def load_ifcb_carbon(data_dir: Path, data_version: str = "fill") -> tuple[pd.DataFrame, list[str]]:
    """Load IFCB carbon biomass and identify taxon biomass columns j."""
    if data_version not in {"clean", "fill"}:
        raise ValueError("data_version must be 'clean' or 'fill'.")

    carbon_path = data_dir / f"ifcb_carbon_{data_version}.csv"
    taxonomy_path = data_dir / f"ifcb_taxonomy_{data_version}.csv" if data_version == "fill" else data_dir / "ifcb_taxonomy.csv"
    if not carbon_path.exists():
        raise FileNotFoundError(
            f"Missing {carbon_path}. Run `ifcb-fill-missing {data_dir}` first, or use --data-version clean."
        )
    if not taxonomy_path.exists():
        raise FileNotFoundError(f"Missing taxonomy file: {taxonomy_path}")

    df = pd.read_csv(carbon_path, low_memory=False)
    taxonomy = pd.read_csv(taxonomy_path)
    taxon_name_col = "Label" if "Label" in taxonomy.columns else "Annotations"
    taxa_cols = list(dict.fromkeys(col for col in taxonomy[taxon_name_col].dropna().astype(str) if col in df.columns))

    # Convert timestamps and taxon biomass to numerical forms used in X[t, i, j].
    df["sample_time"] = pd.to_datetime(df["sample_time"], utc=True, errors="coerce")
    df[taxa_cols] = df[taxa_cols].apply(pd.to_numeric, errors="coerce")

    # Count observations within each cruise x cast; n_obs is a sampling-coverage score.
    df = df.loc[df["year"] != 2026].copy()
    df["n_obs"] = df.groupby(["cruise", "cast"])["pid"].transform("size")

    # JAS 2020 lacks standard L8; use the northernmost L9 underway row as L8.
    stand_in = df.loc[
        (df["year"] == 2020)
        & (df["season"] == "JAS")
        & (df["nearest_station"] == "L9")
        & (df["sample_type"] == "underway")
    ].copy()
    if not stand_in.empty:
        stand_in = stand_in.sort_values("latitude", ascending=False).head(1)
        stand_in["nearest_station"] = "L8"
        stand_in["cast"] = "L8"
        stand_in["sample_type"] = "cast_from_udw"
        df = pd.concat([df, stand_in], ignore_index=True)

    return df, taxa_cols


def select_season_metacommunity(df: pd.DataFrame, season: str) -> pd.DataFrame:
    """Select one balanced surface metacommunity snapshot per year."""
    ds = df.loc[
        df["nearest_station"].isin(STATION_LIST)
        & df["sample_type"].isin(["cast", "cast_from_udw"])
        & (df["season"] == season)
    ].copy()

    # For each year x station x cast, retain the shallowest sample: min(depth).
    ds = ds.sort_values("depth", ascending=True)
    ds = ds.drop_duplicates(["year", "nearest_station", "cast"], keep="first")

    # For each year x station, prefer greatest n_obs, then a known main cruise.
    ds["main_cruise_priority"] = ds["cruise"].isin(MAIN_CRUISE).astype(int)
    ds["station_order"] = pd.Categorical(ds["nearest_station"], categories=STATION_LIST, ordered=True)
    ds = ds.sort_values(
        ["year", "station_order", "n_obs", "main_cruise_priority"],
        ascending=[True, True, False, False],
    )
    ds = ds.drop_duplicates(["year", "nearest_station"], keep="first")

    # Keep complete years only: n_sites(t) = length(STATION_LIST).
    complete_years = (
        ds.groupby("year")["nearest_station"].nunique().loc[lambda count: count == len(STATION_LIST)].index
    )
    ds = ds.loc[ds["year"].isin(complete_years)].copy()
    return ds.drop(columns=["main_cruise_priority", "station_order"])


def build_community_array(df: pd.DataFrame, taxa_cols: list[str]):
    """Build X[t, i, j] where t = year, i = station, and j = taxon."""
    community_wide = df.loc[:, ["nearest_station", "year"] + taxa_cols].rename(
        columns={"nearest_station": "site", "year": "timestep"}
    )
    return make_community_array(community_wide, taxa_cols)


def add_metacommunity_rows(df: pd.DataFrame, taxa_cols: list[str]) -> pd.DataFrame:
    """Add regional rows X_.tj = sum_i X_itj for composition summaries."""
    totals = df.groupby("year")[taxa_cols].sum(min_count=1)
    totals = totals.assign(nearest_station="Metacommunity", sample_time=df.groupby("year")["sample_time"].min())
    totals = totals.reset_index()
    common_cols = [col for col in df.columns if col in totals.columns]
    return pd.concat([df.loc[:, common_cols], totals.loc[:, common_cols]], ignore_index=True)


def dominant_taxa_by_station(df: pd.DataFrame, taxa_cols: list[str], n_per_station: int = 3) -> list[str]:
    """Pool the largest local biomass contributors at each station."""
    totals = df.groupby("nearest_station")[taxa_cols].sum(min_count=1)
    selected: set[str] = set()
    for _, row in totals.iterrows():
        selected.update(row.sort_values(ascending=False).head(n_per_station).index.tolist())
    pooled = totals.loc[:, list(selected)].sum(axis=0).sort_values(ascending=False)
    return pooled.index.tolist()


def prepare_composition_long(df: pd.DataFrame, taxa_cols: list[str], top_species: list[str]) -> pd.DataFrame:
    """Prepare long biomass composition with non-dominant taxa grouped as Others."""
    long = df.melt(
        id_vars=["nearest_station", "year"],
        value_vars=taxa_cols,
        var_name="species",
        value_name="count",
    )
    long["species_grp"] = long["species"].where(long["species"].isin(top_species), "Others")
    out = (
        long.groupby(["nearest_station", "year", "species_grp"], as_index=False)["count"]
        .sum(min_count=1)
        .rename(columns={"nearest_station": "site"})
    )
    return out
