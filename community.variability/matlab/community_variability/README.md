# MATLAB Community Variability

This folder contains a dependency-free MATLAB port of the core
`community.variability` calculations.

The shared data convention is:

```matlab
X(time, site, taxon)
```

You can pass either the numeric array directly or a struct returned by
`make_community_array`, which keeps labels for `timestep`, `site`, and `taxon`.

## Available Functions

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

## Metrics and Notation

The core array convention is `X(time, site, taxon)`. Let `X_{tij}` be biomass
of taxon `j` at time `t` in site `i`.

$$
X_{ti\cdot} = \sum_j X_{tij}, \quad
X_{t\cdot j} = \sum_i X_{tij}, \quad
X_{t\cdot\cdot} = \sum_i \sum_j X_{tij}
$$

$$
\mu_\gamma = \mathrm{mean}_t(X_{t\cdot\cdot}), \quad
\sigma_\gamma = \mathrm{sd}_t(X_{t\cdot\cdot}), \quad
\sigma_i = \mathrm{sd}_t(X_{ti\cdot})
$$

$$
\mu_i = \mathrm{mean}_t(X_{ti\cdot}), \quad
w_i = \frac{\mu_i}{\sum_i \mu_i}, \quad
w_t = \frac{X_{t\cdot\cdot}}{\sum_t X_{t\cdot\cdot}}
$$

$$
CV_\alpha^2 =
\left(
\frac{\sum_i \sigma_i}{\mu_\gamma}
\right)^2,
\quad
CV_\gamma^2 =
\left(
\frac{\sigma_\gamma}{\mu_\gamma}
\right)^2,
\quad
\phi =
\frac{CV_\gamma^2}{CV_\alpha^2}
$$

Compositional variability uses Hellinger composition:

$$
z_{tij} = \sqrt{\frac{X_{tij}}{X_{ti\cdot}}},
\quad
z_{t\cdot j} =
\sqrt{\frac{X_{t\cdot j}}{X_{t\cdot\cdot}}}
$$

$$
\sigma^2_{ij} = \mathrm{Var}_t(z_{tij}), \quad
\sigma^2_{\gamma j} = \mathrm{Var}_t(z_{t\cdot j}), \quad
\sigma^2_{tj} = \mathrm{Var}_i(z_{tij})
$$

$$
BD_\alpha^h =
\sum_i w_i \sum_j \sigma^2_{ij},
\quad
BD_\gamma^h =
\sum_j \sigma^2_{\gamma j},
\quad
BD_\phi^h =
\frac{BD_\gamma^h}{BD_\alpha^h}
$$

$$
BD_\beta^h =
\sum_t w_t \sum_j \sigma^2_{tj}
$$

In MATLAB, `make_community_array()` builds `X(time, site, taxon)`. Metric
functions map directly to the equations with dimension reductions such as
`sum`, `std`, `var`, and `mean` using `"omitnan"`, producing the `\mu`,
`\sigma`, and `\sigma^2` terms.

## Usage

### Compute All Metrics Together

Use this when you want the complete metric table with `CV_alpha`, `CV_gamma`,
`CV_phi`, `BD_alpha`, `BD_gamma`, `BD_phi`, and `BD_beta`.

```matlab
addpath("community.variability/matlab/community_variability")

community = make_community_array(data_wide, taxa_cols);
metrics = calc_metacommunity_metrics(community);
```

`calc_metacommunity_metrics()` returns a table with one row per computed
variable.

### Compute One Metric at a Time

Use this when you want a single computed variable or want to inspect the
spatial component before it is summarized as `BD_beta`.

```matlab
cvAlpha = cv_alpha(community);
cvGamma = cv_gamma(community);
cvPhi = cv_phi(community);

bdAlpha = bd_alpha(community);
bdGamma = bd_gamma(community);
bdPhi = bd_phi(community);
bdBeta = bd_spatial_weighted(community);

spatialByTime = spatial_bd_by_time(community);
```

The individual metric functions return scalar values. `spatial_bd_by_time()`
returns one row per timestep with spatial compositional variability, biomass
weights, and weighted contributions.

## References

Lamy, T. et al. 2021. The dual nature of metacommunity variability. *Oikos*
130: 2078-2092. https://doi.org/10.1111/oik.08517

Git repo: https://github.com/sokole/ltermetacommunities/tree/master/ltmc
