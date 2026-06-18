"""IFCB single-season community-variability workflow in Python."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ifcb_common import (
    add_logging_arguments,
    add_metacommunity_rows,
    build_community_array,
    configure_workflow_logging,
    default_data_dir,
    default_results_dir,
    dominant_taxa_by_station,
    load_ifcb_carbon,
    prepare_composition_long,
    select_season_metacommunity,
)
from community_variability import calc_metacommunity_metrics, wide_metric_table

LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Calculate one seasonal IFCB metacommunity analysis.")
    parser.add_argument("--data-dir", default=default_data_dir(), type=str)
    parser.add_argument("--data-version", choices=["clean", "fill"], default="fill")
    parser.add_argument("--results-dir", default=default_results_dir(), type=str)
    parser.add_argument("--season", default="JAS")
    parser.add_argument("--top-taxa-per-station", default=3, type=int)
    add_logging_arguments(parser)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_workflow_logging(args, "ifcb_single_season")
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Starting single-season workflow for %s", args.season)

    # Load taxon biomass columns j and select one comparable surface sample per year x station.
    df, taxa_cols = load_ifcb_carbon(Path(args.data_dir), data_version=args.data_version)
    ds = select_season_metacommunity(df, args.season)
    LOGGER.info("Complete years retained for %s: %s", args.season, ds["year"].nunique())

    # Build X[t, i, j] and calculate alpha-gamma-phi partitions.
    community_array = build_community_array(ds, taxa_cols)
    metrics = calc_metacommunity_metrics(community_array)
    LOGGER.info("Calculated metrics:\n%s", metrics.to_string(index=False))
    estimate_path = results_dir / f"estimate_{args.season}.csv"
    wide_metric_table(metrics).to_csv(estimate_path, index=False)
    LOGGER.info("Saved metric estimates to: %s", estimate_path)

    # Save composition data used by the R plotting workflow.
    ds_with_meta = add_metacommunity_rows(ds, taxa_cols)
    top_species = dominant_taxa_by_station(ds, taxa_cols, n_per_station=args.top_taxa_per_station)
    composition = prepare_composition_long(ds_with_meta, taxa_cols, top_species)
    composition_path = results_dir / f"composition_{args.season}.csv"
    composition.to_csv(composition_path, index=False)
    LOGGER.info("Saved composition data to: %s", composition_path)
    LOGGER.info("Single-season workflow completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
