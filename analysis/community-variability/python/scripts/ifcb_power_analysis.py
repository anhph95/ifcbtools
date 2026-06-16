"""IFCB seasonal bootstrap power analysis in Python."""

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
from community_variability import bootstrap_by_dimension, calc_spatial_bd_by_time


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Bootstrap seasonal IFCB metacommunity metrics.")
    parser.add_argument("--data-dir", default=default_data_dir(), type=str)
    parser.add_argument("--data-version", choices=["clean", "fill"], default="fill")
    parser.add_argument("--results-dir", default=default_results_dir(), type=str)
    parser.add_argument("--seasons", nargs="+", default=SEASONS)
    parser.add_argument("--n-boot", default=1000, type=int)
    parser.add_argument("--seed", default=123, type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load once; each season applies the same balanced-metacommunity selection.
    df, taxa_cols = load_ifcb_carbon(Path(args.data_dir), data_version=args.data_version)

    for season in args.seasons:
        print(f"Processing season: {season}")
        ds = select_season_metacommunity(df, season)
        print(f"Complete years retained: {ds['year'].nunique()}")

        # Resample timesteps t: X*[t, i, j] = X[sampled t, i, j].
        community_array = build_community_array(ds, taxa_cols)
        boot_res = bootstrap_by_dimension(
            community_array,
            margin="timestep",
            n_boot=args.n_boot,
            seed=args.seed,
            baseline_in_boot=True,
        )
        boot_res["boot"].to_csv(results_dir / f"boot_{season}.csv", index=False)
        boot_res["summary"].to_csv(results_dir / f"summary_{season}.csv", index=False)
        calc_spatial_bd_by_time(community_array).to_csv(
            results_dir / f"spatial_variance_{season}.csv",
            index=False,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
