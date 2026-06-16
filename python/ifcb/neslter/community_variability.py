"""Metacommunity variability metrics for IFCB community arrays."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
import warnings

import numpy as np
import pandas as pd


COMMUNITY_METRIC_ORDER = [
    "CV_gamma",
    "CV_alpha",
    "CV_phi",
    "BD_gamma",
    "BD_alpha",
    "BD_phi",
    "BD_beta",
]


@dataclass(frozen=True)
class CommunityArray:
    """Numerical community array with labels for X[time, site, taxon]."""

    values: np.ndarray
    timestep: list[str]
    site: list[str]
    taxon: list[str]

    @property
    def dimnames(self) -> dict[str, list[str]]:
        return {"timestep": self.timestep, "site": self.site, "taxon": self.taxon}


def _as_array(X: CommunityArray | np.ndarray) -> np.ndarray:
    values = X.values if isinstance(X, CommunityArray) else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")
    return values


def _nanvar_sample(values: np.ndarray, axis: int) -> np.ndarray:
    """Match R stats::var(..., na.rm = TRUE), which uses ddof = 1."""
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanvar(values, axis=axis, ddof=1)


def _nansd_sample(values: np.ndarray, axis: int) -> np.ndarray:
    """Match R stats::sd(..., na.rm = TRUE), which uses ddof = 1."""
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanstd(values, axis=axis, ddof=1)


def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    """Divide and replace non-finite ratios with zero biomass proportion."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.divide(numerator, denominator)
    return np.where(np.isfinite(ratio), ratio, 0.0)


def _safe_scalar_divide(numerator: float, denominator: float) -> float:
    """Match R scalar division behavior for zero denominators."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.divide(numerator, denominator))


def _resolve_margin(X: CommunityArray | np.ndarray, margin: str | int) -> int:
    if isinstance(margin, str):
        names = ["timestep", "site", "taxon"]
        if margin not in names:
            raise ValueError(f"margin must be one of {names}, not {margin!r}.")
        return names.index(margin)
    margin_idx = int(margin)
    if margin_idx not in (0, 1, 2):
        raise ValueError("Integer margin must be 0, 1, or 2.")
    return margin_idx


def _slice_array(X: CommunityArray | np.ndarray, indices: list[object]) -> CommunityArray | np.ndarray:
    values = _as_array(X)[tuple(indices)]
    if not isinstance(X, CommunityArray):
        return values
    labels = X.dimnames
    sliced_labels: dict[str, list[str]] = {}
    for axis, name in enumerate(["timestep", "site", "taxon"]):
        axis_index = indices[axis]
        if isinstance(axis_index, np.ndarray):
            sliced_labels[name] = [labels[name][int(i)] for i in axis_index]
        elif isinstance(axis_index, list):
            sliced_labels[name] = [labels[name][int(i)] for i in axis_index]
        else:
            sliced_labels[name] = list(labels[name])
    return CommunityArray(
        values=values,
        timestep=sliced_labels["timestep"],
        site=sliced_labels["site"],
        taxon=sliced_labels["taxon"],
    )


def make_community_array(
    data_wide: pd.DataFrame,
    taxa_cols: Sequence[str],
    time_step_col_name: str = "timestep",
    site_id_col_name: str = "site",
) -> CommunityArray:
    """Build X[t, i, j] from one wide biomass row per time x site sample."""
    time_ids = sorted(data_wide[time_step_col_name].dropna().unique().tolist())
    site_ids = sorted(data_wide[site_id_col_name].dropna().unique().tolist())
    taxa_cols = list(taxa_cols)

    X = np.zeros((len(time_ids), len(site_ids), len(taxa_cols)), dtype=float)
    time_index = {value: idx for idx, value in enumerate(time_ids)}
    site_index = {value: idx for idx, value in enumerate(site_ids)}

    biomass = data_wide.loc[:, taxa_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float).copy()
    biomass[~np.isfinite(biomass)] = 0.0

    # Fill X[t, i, j] row by row; duplicate time x site rows are summed.
    for row_idx, row in data_wide.reset_index(drop=True).iterrows():
        t = time_index[row[time_step_col_name]]
        i = site_index[row[site_id_col_name]]
        X[t, i, :] += biomass[row_idx, :]

    return CommunityArray(
        values=X,
        timestep=[str(value) for value in time_ids],
        site=[str(value) for value in site_ids],
        taxon=[str(value) for value in taxa_cols],
    )


def cv_gamma(X: CommunityArray | np.ndarray) -> float:
    """Regional aggregate variability: CV_gamma^2 = (sd_t(X_.t.) / mean_t(X_.t.))^2."""
    values = _as_array(X)
    total_metacommunity_biomass = np.nansum(values, axis=(1, 2))
    mu_tt = np.nanmean(total_metacommunity_biomass)
    sigma_tt = _nansd_sample(total_metacommunity_biomass, axis=0)
    return float((sigma_tt / mu_tt) ** 2)


def cv_alpha(X: CommunityArray | np.ndarray) -> float:
    """Local aggregate variability: CV_alpha^2 = (sum_i sd_t(X_it.) / mean_t(X_.t.))^2."""
    values = _as_array(X)
    site_biomass_by_time = np.nansum(values, axis=2)
    total_metacommunity_biomass = np.nansum(site_biomass_by_time, axis=1)
    mu_tt = np.nanmean(total_metacommunity_biomass)
    site_sd = _nansd_sample(site_biomass_by_time, axis=0)
    return float((np.nansum(site_sd) / mu_tt) ** 2)


def cv_phi(X: CommunityArray | np.ndarray) -> float:
    """Spatial aggregate synchrony: phi = CV_gamma^2 / CV_alpha^2."""
    return _safe_scalar_divide(cv_gamma(X), cv_alpha(X))


def bd_gamma(X: CommunityArray | np.ndarray) -> float:
    """Regional compositional variability: BD_gamma^h = sum_j Var_t(z_.tj)."""
    values = _as_array(X)
    regional_taxon_biomass = np.nansum(values, axis=1)
    total_metacommunity_biomass = np.nansum(regional_taxon_biomass, axis=1)

    # Hellinger regional composition: z_.tj = sqrt(X_.tj / X_.t.).
    regional_relative_biomass = _safe_ratio(
        regional_taxon_biomass,
        total_metacommunity_biomass[:, None],
    )
    regional_hellinger = np.sqrt(regional_relative_biomass)
    taxon_temporal_var = _nanvar_sample(regional_hellinger, axis=0)
    return float(np.nansum(taxon_temporal_var))


def bd_alpha(X: CommunityArray | np.ndarray) -> float:
    """Biomass-weighted local compositional variability: BD_alpha^h = sum_i w_i BD_i^h."""
    values = _as_array(X)
    site_biomass_by_time = np.nansum(values, axis=2)

    # Hellinger local composition: z_itj = sqrt(X_itj / X_it.).
    site_relative_biomass = _safe_ratio(values, site_biomass_by_time[:, :, None])
    site_hellinger = np.sqrt(site_relative_biomass)

    # Site-level variability: BD_i^h = sum_j Var_t(z_itj).
    taxon_temporal_var = _nanvar_sample(site_hellinger, axis=0)
    site_bd = np.nansum(taxon_temporal_var, axis=1)

    # Biomass weights: w_i = mean_t(X_it.) / sum_i mean_t(X_it.).
    mean_site_biomass = np.nanmean(site_biomass_by_time, axis=0)
    site_weights = mean_site_biomass / np.nansum(mean_site_biomass)
    return float(np.nansum(site_bd * site_weights))


def bd_phi(X: CommunityArray | np.ndarray) -> float:
    """Spatial compositional synchrony: BD_phi^h = BD_gamma^h / BD_alpha^h."""
    return _safe_scalar_divide(bd_gamma(X), bd_alpha(X))


def spatial_bd_by_time(X: CommunityArray | np.ndarray) -> pd.DataFrame:
    """Return per-timestep spatial compositional variability among sites."""
    values = _as_array(X)
    site_biomass_by_time = np.nansum(values, axis=2)

    # Hellinger local composition: z_itj = sqrt(X_itj / X_it.).
    site_relative_biomass = _safe_ratio(values, site_biomass_by_time[:, :, None])
    site_hellinger = np.sqrt(site_relative_biomass)

    # Spatial variability per timestep: BD_t^h = sum_j Var_i(z_itj).
    taxon_spatial_var = _nanvar_sample(site_hellinger, axis=1)
    bd = np.nansum(taxon_spatial_var, axis=1)

    total_metacommunity_biomass = np.nansum(site_biomass_by_time, axis=1)
    weights = total_metacommunity_biomass / np.nansum(total_metacommunity_biomass)
    timestep = X.timestep if isinstance(X, CommunityArray) else [str(i) for i in range(values.shape[0])]

    return pd.DataFrame(
        {
            "timestep": timestep,
            "BD": bd,
            "total_metacommunity_biomass": total_metacommunity_biomass,
            "weights": weights,
            "BD_x_wt": bd * weights,
        }
    )


def bd_spatial_weighted(X: CommunityArray | np.ndarray) -> float:
    """Biomass-weighted mean spatial compositional variability through time."""
    out = spatial_bd_by_time(X)
    return float(out["BD_x_wt"].sum(skipna=True))


def calc_metacommunity_metrics(X: CommunityArray | np.ndarray) -> pd.DataFrame:
    """Calculate alpha-gamma-phi aggregate and compositional partitions."""
    metrics = {
        "CV_gamma": cv_gamma(X),
        "CV_alpha": cv_alpha(X),
        "BD_gamma": bd_gamma(X),
        "BD_alpha": bd_alpha(X),
        "BD_beta": bd_spatial_weighted(X),
    }
    metrics["CV_phi"] = _safe_scalar_divide(metrics["CV_gamma"], metrics["CV_alpha"])
    metrics["BD_phi"] = _safe_scalar_divide(metrics["BD_gamma"], metrics["BD_alpha"])
    return pd.DataFrame(
        {
            "varname": COMMUNITY_METRIC_ORDER,
            "estimate": [metrics[name] for name in COMMUNITY_METRIC_ORDER],
        }
    )


def wide_metric_table(metric_table: pd.DataFrame, value_col: str = "estimate") -> pd.DataFrame:
    """Convert long metric output to one wide row."""
    return metric_table.pivot_table(columns="varname", values=value_col, aggfunc="first").reset_index(drop=True)


def calc_spatial_bd_by_time(X: CommunityArray | np.ndarray) -> pd.DataFrame:
    """Alias matching the R helper name."""
    return spatial_bd_by_time(X)


def bootstrap_by_dimension(
    X: CommunityArray | np.ndarray,
    margin: str | int,
    metric_fn: Callable[[CommunityArray | np.ndarray], pd.DataFrame] = calc_metacommunity_metrics,
    n_boot: int = 1000,
    seed: int = 123,
    baseline_in_boot: bool = True,
) -> dict[str, pd.DataFrame]:
    """Bootstrap X by resampling one dimension with replacement."""
    rng = np.random.default_rng(seed)
    margin_idx = _resolve_margin(X, margin)
    n_margin = _as_array(X).shape[margin_idx]
    baseline = metric_fn(X)

    boot_rows = []
    for boot_id in range(1, n_boot + 1):
        sampled_idx = rng.integers(0, n_margin, size=n_margin)
        indices: list[object] = [slice(None), slice(None), slice(None)]
        indices[margin_idx] = sampled_idx
        boot_metrics = wide_metric_table(metric_fn(_slice_array(X, indices)))
        boot_metrics.insert(0, "sample_type", "Bootstrap")
        boot_metrics.insert(0, "boot_id", boot_id)
        boot_rows.append(boot_metrics)

    boot_replicates = pd.concat(boot_rows, ignore_index=True) if boot_rows else pd.DataFrame()
    baseline_wide = wide_metric_table(baseline)
    baseline_wide.insert(0, "sample_type", "Baseline")
    baseline_wide.insert(0, "boot_id", 0)
    boot = pd.concat([baseline_wide, boot_replicates], ignore_index=True) if baseline_in_boot else boot_replicates

    summary = baseline.copy()
    summary["lwr"] = [boot_replicates[name].quantile(0.025) for name in summary["varname"]]
    summary["upr"] = [boot_replicates[name].quantile(0.975) for name in summary["varname"]]
    return {
        "baseline": baseline,
        "boot": boot,
        "boot_replicates": boot_replicates,
        "summary": summary,
    }


def leave_one_out(
    X: CommunityArray | np.ndarray,
    margin: str | int,
    metric_fn: Callable[[CommunityArray | np.ndarray], pd.DataFrame] = calc_metacommunity_metrics,
) -> pd.DataFrame:
    """Recalculate metrics after removing each timestep, site, or taxon slice."""
    margin_idx = _resolve_margin(X, margin)
    margin_name = ["timestep", "site", "taxon"][margin_idx]
    removed_col = f"{margin_name}_removed"
    labels = X.dimnames[margin_name] if isinstance(X, CommunityArray) else [str(i) for i in range(_as_array(X).shape[margin_idx])]

    baseline = metric_fn(X)
    baseline.insert(0, removed_col, "Baseline")
    rows = [baseline]

    for label_idx, label in enumerate(labels):
        keep_idx = np.array([i for i in range(len(labels)) if i != label_idx], dtype=int)
        indices: list[object] = [slice(None), slice(None), slice(None)]
        indices[margin_idx] = keep_idx
        omitted_metrics = metric_fn(_slice_array(X, indices))
        omitted_metrics.insert(0, removed_col, label)
        rows.append(omitted_metrics)

    return pd.concat(rows, ignore_index=True)


def add_baseline_delta(
    sensitivity_results: pd.DataFrame,
    removed_col: str,
    baseline_label: str = "Baseline",
    group_cols: str | Sequence[str] | None = None,
) -> pd.DataFrame:
    """Add delta = estimate_without_component - estimate_baseline."""
    group_cols = [] if group_cols is None else ([group_cols] if isinstance(group_cols, str) else list(group_cols))
    join_cols = group_cols + ["varname"]
    baseline = (
        sensitivity_results.loc[sensitivity_results[removed_col] == baseline_label, join_cols + ["estimate"]]
        .rename(columns={"estimate": "baseline"})
        .copy()
    )
    out = sensitivity_results.merge(baseline, on=join_cols, how="left")
    out["delta"] = out["estimate"] - out["baseline"]
    out["abs_delta"] = out["delta"].abs()
    return out


# R-style aliases keep mathematical notation recognizable in Python analyses.
CV_gamma = cv_gamma
CV_alpha = cv_alpha
CV_phi = cv_phi
BD_gamma = bd_gamma
BD_alpha = bd_alpha
BD_phi = bd_phi
BD_spatial_weighted = bd_spatial_weighted
