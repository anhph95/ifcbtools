"""Community array construction and metric output helpers."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .variability import bd_alpha, bd_gamma, bd_spatial_weighted, cv_alpha, cv_gamma


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

    # Initialize the biomass array:
    # X[t, i, j] = biomass of taxon j at site i and time t.
    X = np.zeros((len(time_ids), len(site_ids), len(taxa_cols)), dtype=float)
    time_index = {value: idx for idx, value in enumerate(time_ids)}
    site_index = {value: idx for idx, value in enumerate(site_ids)}

    # Extract taxon biomass columns as a numeric matrix.
    # Non-finite values are treated as zero biomass.
    biomass = data_wide.loc[:, taxa_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float).copy()
    biomass[~np.isfinite(biomass)] = 0.0

    # Fill X row by row. Duplicate site x time rows are added so total
    # biomass is preserved rather than silently dropping replicate records.
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


def calc_metacommunity_metrics(X) -> pd.DataFrame:
    """Calculate alpha-gamma-phi aggregate and compositional partitions."""
    # Calculate alpha and gamma components first.
    # CV_gamma and CV_alpha use aggregate biomass; BD_gamma and BD_alpha use
    # Hellinger composition, z = sqrt(relative biomass).
    metrics = {
        "CV_gamma": cv_gamma(X),
        "CV_alpha": cv_alpha(X),
        "BD_gamma": bd_gamma(X),
        "BD_alpha": bd_alpha(X),
        "BD_beta": bd_spatial_weighted(X),
    }

    # Multiplicative partitions:
    # CV_gamma = CV_alpha * CV_phi
    # BD_gamma = BD_alpha * BD_phi
    with np.errstate(divide="ignore", invalid="ignore"):
        metrics["CV_phi"] = float(np.divide(metrics["CV_gamma"], metrics["CV_alpha"]))
        metrics["BD_phi"] = float(np.divide(metrics["BD_gamma"], metrics["BD_alpha"]))

    return pd.DataFrame(
        {
            "varname": COMMUNITY_METRIC_ORDER,
            "estimate": [metrics[name] for name in COMMUNITY_METRIC_ORDER],
        }
    )


def wide_metric_table(metric_table: pd.DataFrame, value_col: str = "estimate") -> pd.DataFrame:
    """Convert long metric output to one wide row."""
    return metric_table.pivot_table(columns="varname", values=value_col, aggfunc="first").reset_index(drop=True)


def calc_spatial_bd_by_time(X) -> pd.DataFrame:
    """Return time-resolved spatial compositional variability."""
    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    site_biomass_by_time = np.nansum(values, axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        site_relative_biomass = values / site_biomass_by_time[:, :, None]
    site_relative_biomass = np.where(np.isfinite(site_relative_biomass), site_relative_biomass, 0.0)
    site_hellinger = np.sqrt(site_relative_biomass)

    # Spatial variability per timestep:
    # BD_t^h = sum_j Var_i(z_itj).
    taxon_spatial_var = np.nanvar(site_hellinger, axis=1, ddof=1)
    bd = np.nansum(taxon_spatial_var, axis=1)

    total_metacommunity_biomass = np.nansum(site_biomass_by_time, axis=1)
    weights = total_metacommunity_biomass / np.nansum(total_metacommunity_biomass)
    timestep = X.timestep if hasattr(X, "timestep") else [str(i) for i in range(values.shape[0])]

    return pd.DataFrame(
        {
            "timestep": timestep,
            "BD": bd,
            "total_metacommunity_biomass": total_metacommunity_biomass,
            "weights": weights,
            "BD_x_wt": bd * weights,
        }
    )
