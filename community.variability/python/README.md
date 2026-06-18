# Python Community Variability

This installable Python package contains the dependency-light
`community_variability` implementation of the core metacommunity variability
metrics.

## Install

From GitHub without manually cloning the repository:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=community.variability/python"
```

From a local checkout:

```bash
pip install -e community.variability/python
```

## Metrics and Notation

The core array convention is:

```text
X[time, site, taxon]
```

Let `X_{tij}` be biomass of taxon `j` at time `t` in site `i`.

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

In Python, `make_community_array()` builds `X[time, site, taxon]`. Metric
functions map directly to the equations with NumPy reductions: `np.nansum()`
collapses selected dimensions, division forms relative biomass, and
`np.nanstd()` or `np.nanvar()` produces the `\sigma` and `\sigma^2` terms.

## Usage

### Compute All Metrics Together

Use this when you want the complete metric table with `CV_alpha`, `CV_gamma`,
`CV_phi`, `BD_alpha`, `BD_gamma`, `BD_phi`, and `BD_beta`.

```python
from community_variability import make_community_array, calc_metacommunity_metrics

community = make_community_array(data_wide, taxa_cols)
metrics = calc_metacommunity_metrics(community)
```

`calc_metacommunity_metrics()` returns a long table with `varname` and
`estimate` columns.

### Compute One Metric at a Time

Use this when you want a single computed variable or want to inspect the
spatial component before it is summarized as `BD_beta`.

```python
from community_variability import (
    bd_alpha,
    bd_gamma,
    bd_phi,
    bd_spatial_weighted,
    calc_spatial_bd_by_time,
    cv_alpha,
    cv_gamma,
    cv_phi,
)

cv_alpha_value = cv_alpha(community)
cv_gamma_value = cv_gamma(community)
cv_phi_value = cv_phi(community)

bd_alpha_value = bd_alpha(community)
bd_gamma_value = bd_gamma(community)
bd_phi_value = bd_phi(community)
bd_beta_value = bd_spatial_weighted(community)

spatial_by_time = calc_spatial_bd_by_time(community)
```

The individual metric functions return scalar values. `calc_spatial_bd_by_time()`
returns one row per timestep with spatial compositional variability, biomass
weights, and weighted contributions.

## Module Layout

The Python modules mirror the R package source files:

```text
community_variability/
|-- variability.py        # CV/BD metric equations
|-- community_metrics.py  # array construction and metric tables
|-- bootstrap.py          # bootstrap resampling
|-- leave_one_out.py      # sensitivity analysis
`-- core.py               # compatibility re-export
```

## References

Lamy, T. et al. 2021. The dual nature of metacommunity variability. *Oikos*
130: 2078-2092. https://doi.org/10.1111/oik.08517

Git repo: https://github.com/sokole/ltermetacommunities/tree/master/ltmc
