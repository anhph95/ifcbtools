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

Aggregate variability:

- `CV_gamma`: regional aggregate variability
- `CV_alpha`: local aggregate variability
- `CV_phi`: spatial aggregate synchrony

Compositional variability:

- `BD_gamma`: regional compositional variability
- `BD_alpha`: biomass-weighted local compositional variability
- `BD_phi`: spatial compositional synchrony
- `BD_beta`: biomass-weighted spatial compositional variability through time

## Main Helpers

```r
make_community_array(data_wide, taxa_cols)
calc_metacommunity_metrics(X)
bootstrap_by_dimension(X, margin = "timestep")
leave_one_out(X, margin = "taxon")
```
