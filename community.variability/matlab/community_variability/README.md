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

The MATLAB implementation uses the same notation as the top-level community
variability documentation:

$$
X_{tij} = \text{biomass of taxon } j \text{ at time } t \text{ in site } i
$$

$$
X_{ti\cdot} = \sum_j X_{tij}, \quad
X_{t\cdot j} = \sum_i X_{tij}, \quad
X_{t\cdot\cdot} = \sum_i \sum_j X_{tij}
$$

How this works in code:
`make_community_array()` builds `X(time, site, taxon)`. Metric functions use
MATLAB dimension reductions such as `sum`, `std`, `var`, and `mean` with
`"omitnan"` so the code follows the equations directly.

Aggregate variability uses total biomass:

$$
CV_\gamma^2 =
\left(
\frac{\mathrm{sd}_t(X_{t\cdot\cdot})}
{\mathrm{mean}_t(X_{t\cdot\cdot})}
\right)^2
$$

$$
CV_\alpha^2 =
\left(
\frac{\sum_i \mathrm{sd}_t(X_{ti\cdot})}
{\mathrm{mean}_t(X_{t\cdot\cdot})}
\right)^2,
\quad
\phi = \frac{CV_\gamma^2}{CV_\alpha^2}
$$

Compositional variability uses Hellinger composition:

$$
z_{tij} = \sqrt{\frac{X_{tij}}{X_{ti\cdot}}},
\quad
z_{t\cdot j} =
\sqrt{\frac{X_{t\cdot j}}{X_{t\cdot\cdot}}}
$$

$$
BD_\gamma^h =
\sum_j \mathrm{Var}_t(z_{t\cdot j})
$$

$$
BD_\alpha^h =
\sum_i
\left(
\frac{\mathrm{mean}_t(X_{ti\cdot})}
{\sum_i \mathrm{mean}_t(X_{ti\cdot})}
\right)
\sum_j \mathrm{Var}_t(z_{tij})
$$

$$
BD_\phi^h =
\frac{BD_\gamma^h}{BD_\alpha^h}
$$

$$
BD_\beta^h =
\sum_t
\left(
\frac{X_{t\cdot\cdot}}
{\sum_t X_{t\cdot\cdot}}
\right)
\sum_j \mathrm{Var}_i(z_{tij})
$$

## References

Lamy, T. et al. 2021. The dual nature of metacommunity variability. *Oikos*
130: 2078-2092. https://doi.org/10.1111/oik.08517

Git repo: https://github.com/sokole/ltermetacommunities/tree/master/ltmc
