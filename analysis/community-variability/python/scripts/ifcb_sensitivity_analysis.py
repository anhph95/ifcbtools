"""IFCB seasonal leave-one-out sensitivity analysis in Python."""

from __future__ import annotations

import argparse
from pathlib import Path

from ifcb_common import (
    SEASONS,
    build_community_array,
    default_data_dir,
    default_results_dir,
    load_ifcb_carbon,
    select_season_metacommunity,
)
from community_variability import leave_one_out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run leave-one-year and leave-one-taxon IFCB sensitivity analyses.")
    parser.add_argument("--data-dir", default=default_data_dir(), type=str)
    parser.add_argument("--data-version", choices=["clean", "fill"], default="fill")
    parser.add_argument("--results-dir", default=default_results_dir(), type=str)
    parser.add_argument("--seasons", nargs="+", default=SEASONS)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load once; each season repeats the same filtering and completeness rules.
    df, taxa_cols = load_ifcb_carbon(Path(args.data_dir), data_version=args.data_version)

    for season in args.seasons:
        print(f"Processing season: {season}")
        ds = select_season_metacommunity(df, season)
        print(f"Complete years retained: {ds['year'].nunique()}")
        community_array = build_community_array(ds, taxa_cols)

        # Leave-one-year removes X[t, , ]; deltas can be added downstream.
        year_out = leave_one_out(community_array, margin="timestep")
        year_out = year_out.rename(columns={"timestep_removed": "year_removed"})
        year_out.insert(0, "season", season)
        year_out.to_csv(results_dir / f"leave_one_year_out_{season}.csv", index=False)

        # Leave-one-taxon removes X[, , j] and recalculates all metrics.
        taxon_out = leave_one_out(community_array, margin="taxon")
        taxon_out.insert(0, "season", season)
        taxon_out.to_csv(results_dir / f"leave_one_taxon_out_{season}.csv", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
