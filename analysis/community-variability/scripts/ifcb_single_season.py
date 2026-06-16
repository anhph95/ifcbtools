"""IFCB single-season community-variability workflow in Python."""

from __future__ import annotations

import argparse
from pathlib import Path

from ifcb_common import (
    add_metacommunity_rows,
    build_community_array,
    default_data_dir,
    default_results_dir,
    dominant_taxa_by_station,
    load_ifcb_carbon,
    prepare_composition_long,
    select_season_metacommunity,
)
from ifcb.neslter.community_variability import calc_metacommunity_metrics, wide_metric_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Calculate one seasonal IFCB metacommunity analysis.")
    parser.add_argument("--data-dir", default=default_data_dir(), type=str)
    parser.add_argument("--data-version", choices=["clean", "filled"], default="filled")
    parser.add_argument("--results-dir", default=default_results_dir(), type=str)
    parser.add_argument("--season", default="JAS")
    parser.add_argument("--top-taxa-per-station", default=3, type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load taxon biomass columns j and select one comparable surface sample per year x station.
    df, taxa_cols = load_ifcb_carbon(Path(args.data_dir), data_version=args.data_version)
    ds = select_season_metacommunity(df, args.season)
    print(f"Complete years retained for {args.season}: {ds['year'].nunique()}")

    # Build X[t, i, j] and calculate alpha-gamma-phi partitions.
    community_array = build_community_array(ds, taxa_cols)
    metrics = calc_metacommunity_metrics(community_array)
    print(metrics.to_string(index=False))
    wide_metric_table(metrics).to_csv(results_dir / f"estimate_{args.season}.csv", index=False)

    # Save composition data used by the R plotting workflow.
    ds_with_meta = add_metacommunity_rows(ds, taxa_cols)
    top_species = dominant_taxa_by_station(ds, taxa_cols, n_per_station=args.top_taxa_per_station)
    composition = prepare_composition_long(ds_with_meta, taxa_cols, top_species)
    composition.to_csv(results_dir / f"composition_{args.season}.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
