# Community Variability

This folder contains the metacommunity variability code as separate MATLAB, R,
and Python implementations.

The shared ecological array convention is:

```text
X[time, site, taxon]
```

where `time` is repeated surveys, `site` is local communities, and `taxon` is
the IFCB taxon/class biomass dimension.

## Mathematical Notation

### Community Array

Ecological meaning:
`X_{tij}` is biomass of taxon `j` at time `t` in local community or site `i`.

$$
X_{tij}
$$

How this works in code:
`make_community_array()` converts a wide table into `X[time, site, taxon]`.
In MATLAB the same array is indexed as `X(time, site, taxon)`.

### Local Biomass

Ecological meaning:
`X_{ti\cdot}` is total biomass at time `t` in local community `i`, summed over
all taxa.

$$
X_{ti\cdot} = \sum_j X_{tij}
$$

How this works in code:
The metric functions sum the taxon dimension of `X`. In R this is
`apply(X, c(1, 2), sum, na.rm = TRUE)`, in Python this is
`np.nansum(X, axis=2)`, and in MATLAB this is `sum(X, 3, "omitnan")`.

### Regional Taxon Biomass

Ecological meaning:
`X_{t\cdot j}` is regional biomass of taxon `j` at time `t`, summed across all
local communities.

$$
X_{t\cdot j} = \sum_i X_{tij}
$$

How this works in code:
`BD_gamma()` collapses the site dimension before calculating regional
composition through time.

### Total Metacommunity Biomass

Ecological meaning:
`X_{t\cdot\cdot}` is total metacommunity biomass at time `t`, summed over all
local communities and taxa.

$$
X_{t\cdot\cdot} = \sum_i \sum_j X_{tij}
$$

How this works in code:
`CV_gamma()` sums both the site and taxon dimensions to create one regional
biomass time series.

### Local Relative Biomass

Ecological meaning:
`p_{tij}` is the within-site relative biomass of taxon `j`; it removes the
effect of total biomass in site `i` at time `t`.

$$
p_{tij} = \frac{X_{tij}}{X_{ti\cdot}}
$$

How this works in code:
The compositional metrics divide each site-by-time taxon vector by its local
total biomass. Non-finite values from zero-biomass samples are set to zero.

### Local Hellinger Composition

Ecological meaning:
`z_{tij}` is the Hellinger-transformed local composition. It describes
composition without using absolute biomass.

$$
z_{tij} = \sqrt{p_{tij}}
$$

How this works in code:
`BD_alpha()`, `BD_phi()`, `spatial_bd_by_time()`, and
`BD_spatial_weighted()` take the square root of local relative biomass before
calculating variances.

### Regional Hellinger Composition

Ecological meaning:
`z_{t\cdot j}` is the Hellinger-transformed regional composition of taxon `j`
at time `t`.

$$
z_{t\cdot j} =
\sqrt{\frac{X_{t\cdot j}}{X_{t\cdot\cdot}}}
$$

How this works in code:
`BD_gamma()` first sums each taxon across sites, divides by total
metacommunity biomass at each time, then applies the square-root
transformation.

## Metrics

### Regional Aggregate Variability: `CV_gamma`

Ecological meaning:
Regional aggregate variability is temporal variability in total
metacommunity biomass.

$$
CV_\gamma^2 =
\left(
\frac{\mathrm{sd}_t(X_{t\cdot\cdot})}
{\mathrm{mean}_t(X_{t\cdot\cdot})}
\right)^2
$$

How this works in code:
`CV_gamma()` returns the squared temporal coefficient of variation after
summing `X` over sites and taxa.

### Local Aggregate Variability: `CV_alpha`

Ecological meaning:
Local aggregate variability sums the temporal fluctuations of total biomass
within each local community, then scales that sum by mean regional biomass.

$$
CV_\alpha^2 =
\left(
\frac{\sum_i \mathrm{sd}_t(X_{ti\cdot})}
{\mathrm{mean}_t(X_{t\cdot\cdot})}
\right)^2
$$

How this works in code:
`CV_alpha()` sums taxa within each site, calculates one temporal standard
deviation per site, sums those standard deviations, and divides by mean total
metacommunity biomass.

### Spatial Aggregate Synchrony: `CV_phi`

Ecological meaning:
Spatial aggregate synchrony measures how strongly local-community biomass
fluctuations reinforce regional biomass variability.

$$
CV_\gamma^2 = CV_\alpha^2 \phi
$$

$$
\phi = \frac{CV_\gamma^2}{CV_\alpha^2}
$$

How this works in code:
`CV_phi()` or `calc_metacommunity_metrics()` divides `CV_gamma` by `CV_alpha`.
Both values are already squared coefficients of variation.

### Regional Compositional Variability: `BD_gamma`

Ecological meaning:
Regional compositional variability measures temporal change in metacommunity
composition after all local communities are pooled.

$$
BD_\gamma^h =
\sum_j \mathrm{Var}_t(z_{t\cdot j})
$$

How this works in code:
`BD_gamma()` calculates regional Hellinger composition, then sums temporal
variance across taxa.

### Local Compositional Variability: `BD_alpha`

Ecological meaning:
Local compositional variability measures temporal change in composition within
each local community, then averages local values with biomass weights.

$$
BD_i^h =
\sum_j \mathrm{Var}_t(z_{tij})
$$

$$
w_i =
\frac{\mathrm{mean}_t(X_{ti\cdot})}
{\sum_i \mathrm{mean}_t(X_{ti\cdot})}
$$

$$
BD_\alpha^h =
\sum_i w_i BD_i^h
$$

How this works in code:
`BD_alpha()` calculates Hellinger composition for each site, sums temporal
taxon variances within each site, then applies site biomass weights.

### Spatial Compositional Synchrony: `BD_phi`

Ecological meaning:
Spatial compositional synchrony measures how strongly local compositional
trajectories combine into regional compositional variability.

$$
BD_\gamma^h = BD_\alpha^h BD_\phi^h
$$

$$
BD_\phi^h =
\frac{BD_\gamma^h}{BD_\alpha^h}
$$

How this works in code:
`BD_phi()` or `calc_metacommunity_metrics()` divides `BD_gamma` by
`BD_alpha`.

### Spatial Compositional Variability Through Time: `BD_beta`

Ecological meaning:
`BD_beta` measures local-scale spatial compositional variability. At each time,
it asks how different local communities are from one another in Hellinger
composition, then averages those spatial differences through time.

$$
BD_t^h =
\sum_j \mathrm{Var}_i(z_{tij})
$$

$$
w_t =
\frac{X_{t\cdot\cdot}}
{\sum_t X_{t\cdot\cdot}}
$$

$$
BD_\beta^h =
\sum_t w_t BD_t^h
$$

How this works in code:
`spatial_bd_by_time()` returns `BD_t^h`, total metacommunity biomass, weights,
and weighted contributions for each timestep. `BD_spatial_weighted()` returns
the weighted sum and `calc_metacommunity_metrics()` reports it as `BD_beta`.

## Install

Python from GitHub:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=community.variability/python"
```

Python from a local checkout:

```bash
pip install -e community.variability/python
```

R from GitHub:

```r
install.packages("remotes")
remotes::install_git(
  "https://github.com/anhph95/ifcbtools.git",
  subdir = "community.variability/R/community.variability"
)
```

R from a local checkout:

```bash
R CMD INSTALL community.variability/R/community.variability
```

MATLAB:

```matlab
addpath("community.variability/matlab/community_variability")
```

or run:

```matlab
run("community.variability/matlab/install_community_variability.m")
install_community_variability(true)
```

To fetch only the community variability code from Git:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set community.variability
```

## Folders

```text
community.variability/
|-- python/  # installable community_variability Python package
|-- R/       # installable community.variability R package
`-- matlab/  # dependency-free MATLAB functions
```

## References

Lamy, T. et al. 2021. The dual nature of metacommunity variability. *Oikos*
130: 2078-2092. https://doi.org/10.1111/oik.08517

Git repo: https://github.com/sokole/ltermetacommunities/tree/master/ltmc
