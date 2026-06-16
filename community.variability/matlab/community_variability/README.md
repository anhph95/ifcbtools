# MATLAB Community Variability

This folder contains a dependency-free MATLAB port of the core
`community.variability` calculations.

The shared data convention is:

```matlab
X(time, site, taxon)
```

You can pass either the numeric array directly or a struct returned by
`make_community_array`, which keeps labels for `timestep`, `site`, and `taxon`.

## Core Metrics

```matlab
addpath("community.variability/matlab/community_variability")

metrics = calc_metacommunity_metrics(community);
spatial = spatial_bd_by_time(community);
```

Available metric functions:

```text
cv_gamma.m
cv_alpha.m
cv_phi.m
bd_gamma.m
bd_alpha.m
bd_phi.m
bd_spatial_weighted.m
spatial_bd_by_time.m
```

Higher-level utilities:

```text
bootstrap_by_dimension.m
leave_one_out.m
add_baseline_delta.m
```

## Equations

Aggregate variability uses total biomass:

```text
CV_gamma^2 = [sd_t(X_t..) / mean_t(X_t..)]^2
CV_alpha^2 = [sum_i sd_t(X_ti.) / mean_t(X_t..)]^2
CV_phi     = CV_gamma^2 / CV_alpha^2
```

Compositional variability uses Hellinger composition:

```text
p_tij = X_tij / X_ti.
z_tij = sqrt(p_tij)
```

Then:

```text
BD_gamma^h = sum_j Var_t(z_t.j)
BD_alpha^h = sum_i w_i sum_j Var_t(z_tij)
BD_phi^h   = BD_gamma^h / BD_alpha^h
BD_beta^h  = sum_t W_t sum_j Var_i(z_tij)
```
