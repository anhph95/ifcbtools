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

## Metrics

The R implementation uses the same mathematical notation as the top-level
community variability documentation:

$$
X_{tij} = \text{biomass of taxon } j \text{ at time } t \text{ in site } i
$$

$$
X_{ti\cdot} = \sum_j X_{tij}, \quad
X_{t\cdot j} = \sum_i X_{tij}, \quad
X_{t\cdot\cdot} = \sum_i \sum_j X_{tij}
$$

Code mapping:
`make_community_array()` builds `X[time, site, taxon]`. The metric functions
collapse array dimensions with base R operations such as `apply()`, `rowSums()`,
and `sweep()` so the code follows the equations directly.

Aggregate variability:

- `CV_gamma`: regional aggregate variability
- `CV_alpha`: local aggregate variability
- `CV_phi`: spatial aggregate synchrony

Compositional variability:

- `BD_gamma`: regional compositional variability
- `BD_alpha`: biomass-weighted local compositional variability
- `BD_phi`: spatial compositional synchrony
- `BD_beta`: biomass-weighted spatial compositional variability through time

The core equations are:

$$
CV_\gamma^2 =
\left(
\frac{\operatorname{sd}_t(X_{t\cdot\cdot})}
{\operatorname{mean}_t(X_{t\cdot\cdot})}
\right)^2
$$

$$
CV_\alpha^2 =
\left(
\frac{\sum_i \operatorname{sd}_t(X_{ti\cdot})}
{\operatorname{mean}_t(X_{t\cdot\cdot})}
\right)^2,
\quad
\phi = \frac{CV_\gamma^2}{CV_\alpha^2}
$$

$$
z_{tij} = \sqrt{\frac{X_{tij}}{X_{ti\cdot}}},
\quad
z_{t\cdot j} =
\sqrt{\frac{X_{t\cdot j}}{X_{t\cdot\cdot}}}
$$

$$
BD_\gamma^h =
\sum_j \operatorname{Var}_t(z_{t\cdot j})
$$

$$
BD_\alpha^h =
\sum_i
\left(
\frac{\operatorname{mean}_t(X_{ti\cdot})}
{\sum_i \operatorname{mean}_t(X_{ti\cdot})}
\right)
\sum_j \operatorname{Var}_t(z_{tij})
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
\sum_j \operatorname{Var}_i(z_{tij})
$$

## Main Helpers

```r
make_community_array(data_wide, taxa_cols)
calc_metacommunity_metrics(X)
bootstrap_by_dimension(X, margin = "timestep")
leave_one_out(X, margin = "taxon")
```

## References

Lamy, T. et al. 2021. The dual nature of metacommunity variability. *Oikos*
130: 2078-2092. https://doi.org/10.1111/oik.08517

Git repo: https://github.com/sokole/ltermetacommunities/tree/master/ltmc
