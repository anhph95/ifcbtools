"""Bootstrap resampling for community arrays."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from .community_metrics import CommunityArray, calc_metacommunity_metrics, wide_metric_table


def bootstrap_by_dimension(
    X: CommunityArray | np.ndarray,
    margin: str | int,
    metric_fn: Callable[[CommunityArray | np.ndarray], pd.DataFrame] = calc_metacommunity_metrics,
    n_boot: int = 1000,
    seed: int = 123,
    baseline_in_boot: bool = True,
) -> dict[str, pd.DataFrame]:
    """Bootstrap X by resampling one dimension with replacement."""
    values = X.values if isinstance(X, CommunityArray) else X
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("Community data must have shape X[time, site, taxon].")

    if isinstance(margin, str):
        names = ["timestep", "site", "taxon"]
        if margin not in names:
            raise ValueError(f"margin must be one of {names}, not {margin!r}.")
        margin_idx = names.index(margin)
    else:
        margin_idx = int(margin)
        if margin_idx not in (0, 1, 2):
            raise ValueError("Integer margin must be 0, 1, or 2.")

    rng = np.random.default_rng(seed)
    n_margin = values.shape[margin_idx]
    baseline = metric_fn(X)

    boot_rows = []
    for boot_id in range(1, n_boot + 1):
        # Resample slices along one ecological dimension:
        # timestep t, site i, or taxon j.
        sampled_idx = rng.integers(0, n_margin, size=n_margin)
        indices: list[object] = [slice(None), slice(None), slice(None)]
        indices[margin_idx] = sampled_idx
        boot_values = values[tuple(indices)]

        if isinstance(X, CommunityArray):
            labels = X.dimnames
            sliced_labels: dict[str, list[str]] = {}
            for axis, name in enumerate(["timestep", "site", "taxon"]):
                if axis == margin_idx:
                    sliced_labels[name] = [labels[name][int(i)] for i in sampled_idx]
                else:
                    sliced_labels[name] = list(labels[name])
            boot_X: CommunityArray | np.ndarray = CommunityArray(
                values=boot_values,
                timestep=sliced_labels["timestep"],
                site=sliced_labels["site"],
                taxon=sliced_labels["taxon"],
            )
        else:
            boot_X = boot_values

        boot_metrics = wide_metric_table(metric_fn(boot_X))
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
