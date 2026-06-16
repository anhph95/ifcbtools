"""Sensitivity analyses for community arrays."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd

from .community_metrics import CommunityArray, calc_metacommunity_metrics


def leave_one_out(
    X: CommunityArray | np.ndarray,
    margin: str | int,
    metric_fn: Callable[[CommunityArray | np.ndarray], pd.DataFrame] = calc_metacommunity_metrics,
) -> pd.DataFrame:
    """Recalculate metrics after removing each timestep, site, or taxon slice."""
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

    margin_name = ["timestep", "site", "taxon"][margin_idx]
    removed_col = f"{margin_name}_removed"
    labels = X.dimnames[margin_name] if isinstance(X, CommunityArray) else [str(i) for i in range(values.shape[margin_idx])]

    baseline = metric_fn(X)
    baseline.insert(0, removed_col, "Baseline")
    rows = [baseline]

    for label_idx, label in enumerate(labels):
        # Omit one ecological component, then recompute all metrics.
        # delta = metric_without_component - metric_baseline is added later.
        keep_idx = np.array([i for i in range(len(labels)) if i != label_idx], dtype=int)
        indices: list[object] = [slice(None), slice(None), slice(None)]
        indices[margin_idx] = keep_idx
        omitted_values = values[tuple(indices)]

        if isinstance(X, CommunityArray):
            dimnames = X.dimnames
            sliced_labels: dict[str, list[str]] = {}
            for axis, name in enumerate(["timestep", "site", "taxon"]):
                if axis == margin_idx:
                    sliced_labels[name] = [dimnames[name][int(i)] for i in keep_idx]
                else:
                    sliced_labels[name] = list(dimnames[name])
            omitted_X: CommunityArray | np.ndarray = CommunityArray(
                values=omitted_values,
                timestep=sliced_labels["timestep"],
                site=sliced_labels["site"],
                taxon=sliced_labels["taxon"],
            )
        else:
            omitted_X = omitted_values

        omitted_metrics = metric_fn(omitted_X)
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
