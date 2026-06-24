# Community Variability Analysis

This folder contains analysis workflows for seasonal IFCB metacommunity
variability.

## Inputs

The scripts expect processed Python outputs in:

```text
data/NESLTER_transect/
```

Create those files with the processing package:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=ifcb.process"
ifcb-process data/NESLTER_transect/ifcb_carbon.csv --clean
ifcb-fill-missing data/NESLTER_transect/ifcb_carbon_clean.csv
```

The analysis scripts use the separate community variability packages:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=community.variability/python"
```

```r
remotes::install_git(
  "https://github.com/anhph95/ifcbtools.git",
  subdir = "community.variability/R/community.variability"
)
```

## Fetch Only This Analysis

The analysis folder is workflow code rather than an installable package.
Install the metric package separately, then use Git only to copy the analysis
source:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set analysis/community-variability
```

Use clean data for observed-sample analyses:

```text
ifcb_carbon_clean.csv
ifcb_taxonomy.csv
```

Use fill data for balanced metacommunity analyses:

```text
ifcb_carbon_clean_fill.csv
ifcb_taxonomy.csv
```

## R Workflows

From this directory:

```bash
Rscript R/scripts/ifcb_single_season.R
Rscript R/scripts/ifcb_power_analysis.R
Rscript R/scripts/ifcb_sensitivity_analysis.R
Rscript R/scripts/ifcb_seasonal_comparison.R
```

## Statistical Workflows

These workflows use the metacommunity array:

$$
X_{tij} = \text{biomass of taxon } j \text{ in site } i \text{ at time } t
$$

where each season is analyzed as a balanced annual metacommunity.

### Bootstrap Power and Uncertainty Analysis

Description:
The bootstrap resamples timesteps with replacement. Each selected timestep
keeps its full site-by-taxon community matrix, preserving spatial and
compositional structure within each sampled year.

$$
X^*_{tij} = X_{s_t i j}, \quad s_t \in \{1, \ldots, n_t\}
$$

How it works in code:
`bootstrap_by_dimension(X, margin = "timestep")` recalculates all variability
metrics for each bootstrap replicate, writes replicate-level outputs, and
summarizes 95% intervals from bootstrap quantiles.

Usage:

```bash
Rscript R/scripts/ifcb_power_analysis.R
python analysis/community-variability/python/scripts/ifcb_power_analysis.py
```

```matlab
run("matlab/scripts/ifcb_power_analysis.m")
```

### Leave-One-Year Sensitivity

Description:
Leave-one-year sensitivity removes one timestep and recalculates the metrics.
Large departures from the baseline identify years that strongly influence
seasonal variability estimates.

$$
X^{(-t)} = X[\{1,\ldots,n_t\} \setminus t, i, j]
$$

$$
\Delta_m^{(-t)} = m(X^{(-t)}) - m(X)
$$

How it works in code:
`leave_one_out(X, margin = "timestep")` produces one baseline row and one
metric table for every omitted year. `add_baseline_delta()` adds the change
from baseline.

Usage:

```bash
Rscript R/scripts/ifcb_sensitivity_analysis.R
python analysis/community-variability/python/scripts/ifcb_sensitivity_analysis.py
```

```matlab
run("matlab/scripts/ifcb_sensitivity_analysis.m")
```

### Leave-One-Taxon Sensitivity

Description:
Leave-one-taxon sensitivity removes one taxon and recalculates the metrics.
Large departures from the baseline identify taxa that strongly influence
aggregate variability, compositional variability, or synchrony.

$$
X^{(-j)} = X[t, i, \{1,\ldots,n_j\} \setminus j]
$$

$$
\Delta_m^{(-j)} = m(X^{(-j)}) - m(X)
$$

How it works in code:
`leave_one_out(X, margin = "taxon")` recalculates all metrics after removing
each taxon. The seasonal comparison workflow classifies unusually large
positive or negative taxon effects using empirical 5% and 95% quantiles of
the taxon-removal deltas within each season and metric.

Usage:

```bash
Rscript R/scripts/ifcb_sensitivity_analysis.R
python analysis/community-variability/python/scripts/ifcb_sensitivity_analysis.py
```

```matlab
run("matlab/scripts/ifcb_sensitivity_analysis.m")
```

### Seasonal Metric Comparison

Description:
Seasonal comparison plots combine observed seasonal metric estimates with
bootstrap uncertainty clouds. For each two-metric comparison, a 95% covariance
ellipse summarizes the bootstrap distribution.

$$
(\mathbf{x} - \boldsymbol{\mu})^\top
\mathbf{S}^{-1}
(\mathbf{x} - \boldsymbol{\mu})
=
\chi^2_{0.95,2}
$$

How it works in code:
`ifcb_seasonal_comparison.R` reads bootstrap outputs, estimates the covariance
matrix for each season and metric pair, draws 95% ellipses, and overlays the
observed seasonal estimates.

Usage:

```bash
Rscript R/scripts/ifcb_seasonal_comparison.R
python analysis/community-variability/python/scripts/ifcb_seasonal_comparison.py
```

```matlab
run("matlab/scripts/ifcb_seasonal_comparison.m")
```

### Spatial Compositional Variability Through Time

Description:
This workflow calculates local-scale spatial compositional variability at each
timestep and tracks how among-site compositional differences change through
time.

$$
BD_t^h =
\sum_j \mathrm{Var}_i(z_{tij})
$$

$$
BD_\beta^h =
\sum_t
\left(
\frac{X_{t\cdot\cdot}}
{\sum_t X_{t\cdot\cdot}}
\right)
BD_t^h
$$

How it works in code:
`calc_spatial_bd_by_time()` or `spatial_bd_by_time()` writes one row per
timestep with `BD`, total metacommunity biomass, timestep weight, and weighted
contribution. The seasonal comparison script plots these outputs as time
series.

## Python Workflows

Use the default carbon input:

```bash
python analysis/community-variability/python/scripts/ifcb_single_season.py
python analysis/community-variability/python/scripts/ifcb_power_analysis.py
python analysis/community-variability/python/scripts/ifcb_sensitivity_analysis.py
python analysis/community-variability/python/scripts/ifcb_seasonal_comparison.py
```

Use another carbon input explicitly:

```bash
python analysis/community-variability/python/scripts/ifcb_single_season.py --input-file path/to/carbon.csv
```

## MATLAB Workflows

From MATLAB, run from this directory:

```matlab
run("matlab/install_analysis_path.m")
install_analysis_path(true)
```

```matlab
run("matlab/scripts/ifcb_single_season.m")
run("matlab/scripts/ifcb_power_analysis.m")
run("matlab/scripts/ifcb_sensitivity_analysis.m")
run("matlab/scripts/ifcb_seasonal_comparison.m")
```

## Outputs

Generated files are written to:

```text
analysis/community-variability/results/
```

The results directory is ignored by Git.

## References

Lamy, T. et al. 2021. The dual nature of metacommunity variability. *Oikos*
130: 2078-2092. https://doi.org/10.1111/oik.08517

Git repo: https://github.com/sokole/ltermetacommunities/tree/master/ltmc
