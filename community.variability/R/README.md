# R Package

The R `community.variability` package provides metacommunity variability metrics
for arrays:

```r
X[time, site, taxon]
```

where:

- `time` is repeated surveys `t`
- `site` is local communities `i`
- `taxon` is taxa `j`

## Install

From GitHub without manually cloning the repository:

```r
install.packages("remotes")
remotes::install_git(
  "https://github.com/anhph95/ifcbtools.git",
  subdir = "community.variability/R/community.variability"
)
```

From a local checkout:

```bash
R CMD INSTALL community.variability/R/community.variability
```

or from R:

```r
install.packages("community.variability/R/community.variability", repos = NULL, type = "source")
```

## Metrics and Notation

The R implementation follows the top-level community variability notation.
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

Aggregate variability:

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

In R, `make_community_array()` builds `X[time, site, taxon]`. The metric
functions map directly to the equations with base R array operations:
`apply()` collapses selected dimensions, `sweep()` divides by local or
metacommunity biomass, and `stats::sd()` or `stats::var()` produces the
`\sigma` and `\sigma^2` terms.

Computed variables:

- `CV_alpha`: local aggregate variability
- `CV_gamma`: metacommunity aggregate variability
- `CV_phi`: aggregate synchrony
- `BD_alpha`: biomass-weighted local compositional variability
- `BD_gamma`: metacommunity compositional variability
- `BD_phi`: compositional synchrony
- `BD_beta`: local-scale spatial compositional variability

## Usage

### Compute All Metrics Together

Use this when you want the complete metric table with `CV_alpha`, `CV_gamma`,
`CV_phi`, `BD_alpha`, `BD_gamma`, `BD_phi`, and `BD_beta`.

```r
library(community.variability)

community <- make_community_array(data_wide, taxa_cols)
metrics <- calc_metacommunity_metrics(community)
```

`calc_metacommunity_metrics()` returns a long table with `varname` and
`estimate` columns.

### Compute One Metric at a Time

Use this when you want a single computed variable or want to inspect the
spatial component before it is summarized as `BD_beta`.

```r
cv_alpha_value <- CV_alpha(community)
cv_gamma_value <- CV_gamma(community)
cv_phi_value <- CV_phi(community)

bd_alpha_value <- BD_alpha(community)
bd_gamma_value <- BD_gamma(community)
bd_phi_value <- BD_phi(community)
bd_beta_value <- BD_spatial_weighted(community)

spatial_by_time <- calc_spatial_bd_by_time(community)
```

The individual metric functions return scalar values. `calc_spatial_bd_by_time()`
returns one row per timestep with spatial compositional variability, biomass
weights, and weighted contributions.

## Main Helpers

```r
make_community_array(data_wide, taxa_cols)
calc_metacommunity_metrics(X)
CV_alpha(X)
CV_gamma(X)
CV_phi(X)
BD_alpha(X)
BD_gamma(X)
BD_phi(X)
BD_spatial_weighted(X)
calc_spatial_bd_by_time(X)
bootstrap_by_dimension(X, margin = "timestep")
leave_one_out(X, margin = "taxon")
```

## References

Lamy, T. et al. 2021. The dual nature of metacommunity variability. *Oikos*
130: 2078-2092. https://doi.org/10.1111/oik.08517

Git repo: https://github.com/sokole/ltermetacommunities/tree/master/ltmc
