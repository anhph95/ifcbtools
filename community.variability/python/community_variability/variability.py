"""Metacommunity variability metrics following Lamy et al. 2021.

Input:
    X[time, site, taxon]

Ecological meaning:
    time  = repeated surveys t
    site  = local communities i
    taxon = species/taxa j

Aggregate variability metrics use total biomass. Compositional variability
metrics use Hellinger-transformed relative biomass so that compositional change
is separated from variation in total biomass.
"""

from __future__ import annotations

import warnings

import numpy as np


def cv_gamma(X) -> float:
    """Regional aggregate variability: CV_gamma^2 = (sd_t(X_.t.) / mean_t(X_.t.))^2."""
    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    # Collapse the full metacommunity to regional biomass through time:
    # X_.t. = sum_i sum_j X_itj.
    total_metacommunity_biomass = np.nansum(values, axis=(1, 2))
    mu_tt = np.nanmean(total_metacommunity_biomass)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sigma_tt = np.nanstd(total_metacommunity_biomass, axis=0, ddof=1)

    # CV_gamma^2 is the squared temporal coefficient of variation.
    return float((sigma_tt / mu_tt) ** 2)


def cv_alpha(X) -> float:
    """Local aggregate variability: CV_alpha^2 = (sum_i sd_t(X_it.) / mean_t(X_.t.))^2."""
    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    # Collapse taxa within each local community:
    # X_it. = sum_j X_itj.
    site_biomass_by_time = np.nansum(values, axis=2)
    total_metacommunity_biomass = np.nansum(site_biomass_by_time, axis=1)
    mu_tt = np.nanmean(total_metacommunity_biomass)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        site_sd = np.nanstd(site_biomass_by_time, axis=0, ddof=1)

    # Sum local biomass standard deviations before scaling by regional mean.
    return float((np.nansum(site_sd) / mu_tt) ** 2)


def cv_phi(X) -> float:
    """Spatial aggregate synchrony: phi = CV_gamma^2 / CV_alpha^2."""
    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    total_metacommunity_biomass = np.nansum(values, axis=(1, 2))
    mu_tt = np.nanmean(total_metacommunity_biomass)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sigma_tt = np.nanstd(total_metacommunity_biomass, axis=0, ddof=1)
    cv_g = (sigma_tt / mu_tt) ** 2

    site_biomass_by_time = np.nansum(values, axis=2)
    total_metacommunity_biomass = np.nansum(site_biomass_by_time, axis=1)
    mu_tt = np.nanmean(total_metacommunity_biomass)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        site_sd = np.nanstd(site_biomass_by_time, axis=0, ddof=1)
    cv_a = (np.nansum(site_sd) / mu_tt) ** 2

    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.divide(cv_g, cv_a))


def bd_gamma(X) -> float:
    """Regional compositional variability: BD_gamma^h = sum_j Var_t(z_.tj)."""
    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    # Regional taxon biomass: X_.tj = sum_i X_itj.
    regional_taxon_biomass = np.nansum(values, axis=1)
    total_metacommunity_biomass = np.nansum(regional_taxon_biomass, axis=1)

    # Hellinger regional composition:
    # z_.tj = sqrt(X_.tj / X_.t.).
    with np.errstate(divide="ignore", invalid="ignore"):
        regional_relative_biomass = regional_taxon_biomass / total_metacommunity_biomass[:, None]
    regional_relative_biomass = np.where(np.isfinite(regional_relative_biomass), regional_relative_biomass, 0.0)
    regional_hellinger = np.sqrt(regional_relative_biomass)

    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        taxon_temporal_var = np.nanvar(regional_hellinger, axis=0, ddof=1)
    return float(np.nansum(taxon_temporal_var))


def bd_alpha(X) -> float:
    """Biomass-weighted local compositional variability: BD_alpha^h = sum_i w_i BD_i^h."""
    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    # Local total biomass: X_it. = sum_j X_itj.
    site_biomass_by_time = np.nansum(values, axis=2)

    # Hellinger local composition:
    # z_itj = sqrt(X_itj / X_it.).
    with np.errstate(divide="ignore", invalid="ignore"):
        site_relative_biomass = values / site_biomass_by_time[:, :, None]
    site_relative_biomass = np.where(np.isfinite(site_relative_biomass), site_relative_biomass, 0.0)
    site_hellinger = np.sqrt(site_relative_biomass)

    # Site-level variability: BD_i^h = sum_j Var_t(z_itj).
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        taxon_temporal_var = np.nanvar(site_hellinger, axis=0, ddof=1)
    site_bd = np.nansum(taxon_temporal_var, axis=1)

    # Biomass weights: w_i = mean_t(X_it.) / sum_i mean_t(X_it.).
    mean_site_biomass = np.nanmean(site_biomass_by_time, axis=0)
    site_weights = mean_site_biomass / np.nansum(mean_site_biomass)
    return float(np.nansum(site_bd * site_weights))


def bd_phi(X) -> float:
    """Spatial compositional synchrony: BD_phi^h = BD_gamma^h / BD_alpha^h."""
    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    regional_taxon_biomass = np.nansum(values, axis=1)
    total_metacommunity_biomass = np.nansum(regional_taxon_biomass, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        regional_relative_biomass = regional_taxon_biomass / total_metacommunity_biomass[:, None]
    regional_relative_biomass = np.where(np.isfinite(regional_relative_biomass), regional_relative_biomass, 0.0)
    regional_hellinger = np.sqrt(regional_relative_biomass)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        taxon_temporal_var = np.nanvar(regional_hellinger, axis=0, ddof=1)
    bd_g = np.nansum(taxon_temporal_var)

    site_biomass_by_time = np.nansum(values, axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        site_relative_biomass = values / site_biomass_by_time[:, :, None]
    site_relative_biomass = np.where(np.isfinite(site_relative_biomass), site_relative_biomass, 0.0)
    site_hellinger = np.sqrt(site_relative_biomass)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        taxon_temporal_var = np.nanvar(site_hellinger, axis=0, ddof=1)
    site_bd = np.nansum(taxon_temporal_var, axis=1)
    mean_site_biomass = np.nanmean(site_biomass_by_time, axis=0)
    site_weights = mean_site_biomass / np.nansum(mean_site_biomass)
    bd_a = np.nansum(site_bd * site_weights)

    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.divide(bd_g, bd_a))


def spatial_bd_by_time(X):
    """Return per-timestep spatial compositional variability among sites."""
    import pandas as pd

    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    site_biomass_by_time = np.nansum(values, axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        site_relative_biomass = values / site_biomass_by_time[:, :, None]
    site_relative_biomass = np.where(np.isfinite(site_relative_biomass), site_relative_biomass, 0.0)
    site_hellinger = np.sqrt(site_relative_biomass)

    # Spatial variability per timestep: BD_t^h = sum_j Var_i(z_itj).
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
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


def bd_spatial_weighted(X) -> float:
    """Biomass-weighted mean spatial compositional variability through time."""
    values = X.values if hasattr(X, "values") else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    site_biomass_by_time = np.nansum(values, axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        site_relative_biomass = values / site_biomass_by_time[:, :, None]
    site_relative_biomass = np.where(np.isfinite(site_relative_biomass), site_relative_biomass, 0.0)
    site_hellinger = np.sqrt(site_relative_biomass)

    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        taxon_spatial_var = np.nanvar(site_hellinger, axis=1, ddof=1)
    bd = np.nansum(taxon_spatial_var, axis=1)

    total_metacommunity_biomass = np.nansum(site_biomass_by_time, axis=1)
    weights = total_metacommunity_biomass / np.nansum(total_metacommunity_biomass)
    return float(np.nansum(bd * weights))


# R-style aliases keep mathematical notation recognizable in Python analyses.
CV_gamma = cv_gamma
CV_alpha = cv_alpha
CV_phi = cv_phi
BD_gamma = bd_gamma
BD_alpha = bd_alpha
BD_phi = bd_phi
BD_spatial_weighted = bd_spatial_weighted
