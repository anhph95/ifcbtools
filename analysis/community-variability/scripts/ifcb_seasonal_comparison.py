"""IFCB seasonal comparison summaries in Python."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ifcb_common import SEASONS, default_results_dir
from ifcb.neslter.community_variability import COMMUNITY_METRIC_ORDER, add_baseline_delta


PLOT_VARS_SENSITIVITY = ["BD_gamma", "CV_gamma", "BD_phi", "CV_phi"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare seasonal comparison CSVs from Python analysis outputs.")
    parser.add_argument("--results-dir", default=default_results_dir(), type=str)
    parser.add_argument("--seasons", nargs="+", default=SEASONS)
    parser.add_argument("--top-n-taxa", default=10, type=int)
    return parser.parse_args()


def read_seasonal_csv(results_dir: Path, seasons: list[str], prefix: str) -> pd.DataFrame:
    rows = []
    for season in seasons:
        path = results_dir / f"{prefix}_{season}.csv"
        data = pd.read_csv(path)
        if "season" not in data.columns:
            data.insert(0, "season", season)
        rows.append(data)
    return pd.concat(rows, ignore_index=True)


def make_ellipse(data: pd.DataFrame, seasons: list[str], xvar: str, yvar: str) -> pd.DataFrame:
    """Calculate 95% covariance ellipse points: (x - mu)' S^-1 (x - mu) = chi2_0.95."""
    theta = np.linspace(0.0, 2.0 * np.pi, 200)
    r2 = 5.991464547107979  # qchisq(0.95, df = 2)
    rows = []
    for season in seasons:
        d = data.loc[data["season"] == season, [xvar, yvar]].dropna()
        if len(d) < 3:
            continue
        cov = np.cov(d.to_numpy(dtype=float), rowvar=False, ddof=1)
        if not np.all(np.isfinite(cov)) or np.linalg.det(cov) <= 0:
            continue
        mu = d.mean(axis=0).to_numpy(dtype=float)
        eigval, eigvec = np.linalg.eigh(cov)
        transform = eigvec @ np.diag(np.sqrt(np.maximum(eigval, 0.0) * r2))
        points = mu + np.column_stack([np.cos(theta), np.sin(theta)]) @ transform.T
        rows.append(pd.DataFrame({"season": season, "x": points[:, 0], "y": points[:, 1], "xvar": xvar, "yvar": yvar}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["season", "x", "y", "xvar", "yvar"])


def classify_taxon_effects(loo_taxon_all: pd.DataFrame, top_n_taxa: int) -> pd.DataFrame:
    """Classify taxon removals using seasonal 5th and 95th percentile delta thresholds."""
    thresholds = (
        loo_taxon_all.loc[loo_taxon_all["taxon_removed"] != "Baseline"]
        .groupby(["season", "varname"])["delta"]
        .quantile([0.05, 0.95])
        .unstack()
        .rename(columns={0.05: "lwr", 0.95: "upr"})
        .reset_index()
    )
    data = loo_taxon_all.loc[
        (loo_taxon_all["taxon_removed"] != "Baseline")
        & loo_taxon_all["varname"].isin(PLOT_VARS_SENSITIVITY)
    ].merge(thresholds, on=["season", "varname"], how="left")

    # Positive means removing the taxon increases the metric; negative means it decreases it.
    significant = (data["delta"] < data["lwr"]) | (data["delta"] > data["upr"])
    data["effect_class"] = "Not significant"
    data.loc[significant & (data["delta"] > 0), "effect_class"] = "Positive"
    data.loc[significant & (data["delta"] < 0), "effect_class"] = "Negative"

    return (
        data.sort_values(["season", "varname", "abs_delta"], ascending=[True, True, False])
        .groupby(["season", "varname"], as_index=False)
        .head(top_n_taxa)
        .reset_index(drop=True)
    )


def main() -> int:
    args = parse_args()
    results_dir = Path(args.results_dir)
    seasons = list(args.seasons)

    # Read bootstrap summaries and replicate-level metric tables.
    summary_all = read_seasonal_csv(results_dir, seasons, "summary")
    boot_all = read_seasonal_csv(results_dir, seasons, "boot")
    boot_all = boot_all.loc[boot_all["sample_type"] == "Bootstrap"].copy()

    summary_all["varname"] = pd.Categorical(summary_all["varname"], categories=COMMUNITY_METRIC_ORDER, ordered=True)
    summary_wide = summary_all.pivot(index="season", columns="varname", values=["estimate", "lwr", "upr"])
    summary_wide.columns = [f"{metric}_{value}" for value, metric in summary_wide.columns]
    summary_wide = summary_wide.reset_index()
    summary_wide.to_csv(results_dir / "seasonal_metric_summary_wide.csv", index=False)

    # Leave-one-year and leave-one-taxon deltas: delta = metric_without_component - baseline.
    loo_year_all = read_seasonal_csv(results_dir, seasons, "leave_one_year_out")
    loo_year_all = add_baseline_delta(loo_year_all, removed_col="year_removed", group_cols="season")
    loo_year_all.to_csv(results_dir / "leave_one_year_all_seasons_with_delta.csv", index=False)

    loo_taxon_all = read_seasonal_csv(results_dir, seasons, "leave_one_taxon_out")
    loo_taxon_all = add_baseline_delta(loo_taxon_all, removed_col="taxon_removed", group_cols="season")
    loo_taxon_all.to_csv(results_dir / "leave_one_taxon_all_seasons_with_delta.csv", index=False)
    classify_taxon_effects(loo_taxon_all, args.top_n_taxa).to_csv(
        results_dir / "leave_one_taxon_top_effects_all_seasons.csv",
        index=False,
    )

    # Bootstrap metric-space ellipses used by the R seasonal comparison figures.
    ellipse_pairs = [
        ("BD_alpha", "CV_alpha"),
        ("BD_phi", "CV_phi"),
        ("BD_gamma", "CV_gamma"),
        ("BD_beta", "BD_gamma"),
    ]
    ellipses = [make_ellipse(boot_all, seasons, xvar, yvar) for xvar, yvar in ellipse_pairs]
    pd.concat(ellipses, ignore_index=True).to_csv(results_dir / "bootstrap_metric_ellipses.csv", index=False)

    # Time-resolved spatial compositional variability: BD_t^h = sum_j Var_i(z_itj).
    spavar_all = read_seasonal_csv(results_dir, seasons, "spatial_variance")
    spavar_all["year"] = pd.to_numeric(spavar_all["timestep"], errors="coerce")
    spavar_all["season_index"] = spavar_all["season"].map({season: i for i, season in enumerate(seasons)})
    spavar_all["time_index"] = spavar_all["year"] + spavar_all["season_index"] / len(seasons)
    spavar_all.to_csv(results_dir / "spatial_variance_all_seasons.csv", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
