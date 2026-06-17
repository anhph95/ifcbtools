# Community Variability

This folder contains the metacommunity variability code as separate MATLAB, R,
and Python implementations.

The shared ecological array convention is:

```text
X[time, site, taxon]
```

where `time` is repeated surveys, `site` is local communities, and `taxon` is
the IFCB taxon/class biomass dimension.

## Metrics

The metrics follow Lamy et al. (2021). Let `X_{tij}` be biomass of taxon `j`
at time `t` in local community `i`. Dots indicate summation over an axis, so
`X_{ti\cdot}` is local aggregate biomass and `X_{t\cdot\cdot}` is
metacommunity aggregate biomass.

### `CV_alpha`: Local Aggregate Variability

`CV_alpha` is temporal variability in aggregate biomass at the local-community
scale. It describes how much total biomass fluctuates within sites before
considering whether different sites fluctuate together.

$$
CV_\alpha^2 =
\left(
\frac{\sum_i \mathrm{sd}_t(X_{ti\cdot})}
{\mathrm{mean}_t(X_{t\cdot\cdot})}
\right)^2
$$

In code, aggregate local biomass is represented by summing the taxon dimension
of `X[time, site, taxon]`:

$$
X_{ti\cdot} = \sum_j X_{tij}
$$

The site dimension is retained, so each site has its own biomass time series.
The code calculates a temporal standard deviation for each site, sums those
local standard deviations, and scales the result by mean total metacommunity
biomass.

### `CV_gamma`: Metacommunity Aggregate Variability

`CV_gamma` is temporal variability in aggregate biomass at the metacommunity
scale. It describes how much total biomass fluctuates after all local
communities and taxa are pooled.

$$
CV_\gamma^2 =
\left(
\frac{\mathrm{sd}_t(X_{t\cdot\cdot})}
{\mathrm{mean}_t(X_{t\cdot\cdot})}
\right)^2
$$

In code, total metacommunity biomass is represented by summing both the site
and taxon dimensions:

$$
X_{t\cdot\cdot} = \sum_i \sum_j X_{tij}
$$

This produces one regional biomass time series. `CV_gamma()` returns its
squared temporal coefficient of variation.

### `CV_phi`: Aggregate Synchrony

`CV_phi` is aggregate synchrony. It compares metacommunity aggregate
variability with local aggregate variability to quantify how strongly local
biomass fluctuations combine at the regional scale.

$$
CV_\gamma^2 = CV_\alpha^2 \phi
$$

$$
\phi = \frac{CV_\gamma^2}{CV_\alpha^2}
$$

In code, `CV_phi` is computed as the ratio of `CV_gamma` to `CV_alpha`. Both
stored values are already squared coefficients of variation, so the ratio is
the aggregate synchrony component of the temporal alpha-gamma partition.

### `BD_alpha`: Local Compositional Variability

`BD_alpha` is temporal compositional variability at the local-community scale.
It describes how much taxon composition changes through time within sites,
after removing changes in total biomass.

$$
BD_i^h =
\sum_j \mathrm{Var}_t(z_{tij})
$$

$$
BD_\alpha^h =
\sum_i w_i BD_i^h
$$

In code, each site-time taxon vector is converted to local relative biomass and
then Hellinger composition:

$$
z_{tij} = \sqrt{\frac{X_{tij}}{X_{ti\cdot}}}
$$

For each site, the code calculates temporal variance in `z_{tij}` for every
taxon and sums those taxon variances to obtain `BD_i^h`. Site values are then
averaged using biomass weights,

$$
w_i =
\frac{\mathrm{mean}_t(X_{ti\cdot})}
{\sum_i \mathrm{mean}_t(X_{ti\cdot})}
$$

so sites with more biomass contribute more to `BD_alpha`.

### `BD_gamma`: Metacommunity Compositional Variability

`BD_gamma` is temporal compositional variability at the metacommunity scale.
It describes how regional taxon composition changes through time after pooling
local communities.

$$
BD_\gamma^h =
\sum_j \mathrm{Var}_t(z_{t\cdot j})
$$

In code, each taxon is first summed across sites to produce regional taxon
biomass:

$$
X_{t\cdot j} = \sum_i X_{tij}
$$

The regional taxon vector is then divided by total metacommunity biomass and
Hellinger-transformed:

$$
z_{t\cdot j} =
\sqrt{\frac{X_{t\cdot j}}{X_{t\cdot\cdot}}}
$$

`BD_gamma()` sums the temporal variances of these regional Hellinger taxon
trajectories.

### `BD_phi`: Compositional Synchrony

`BD_phi` is compositional synchrony. It compares metacommunity compositional
variability with local compositional variability to quantify how strongly local
composition trajectories align through time.

$$
BD_\gamma^h = BD_\alpha^h BD_\phi^h
$$

$$
BD_\phi^h =
\frac{BD_\gamma^h}{BD_\alpha^h}
$$

In code, `BD_phi` is computed as the ratio of `BD_gamma` to `BD_alpha`. The
ratio is the compositional synchrony component of the temporal alpha-gamma
partition.

### `BD_beta`: Local-Scale Spatial Compositional Variability

`BD_beta` is compositional variability at the local spatial scale. It measures
how different local communities are from one another in composition at each
time, then summarizes that spatial variability through time.

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

In code, the same local Hellinger composition `z_{tij}` used for `BD_alpha` is
used here, but the variance is taken across sites instead of across time. For
each timestep, `spatial_bd_by_time()` calculates among-site variance for each
taxon and sums those variances to obtain `BD_t^h`. `BD_spatial_weighted()`
weights each `BD_t^h` by total metacommunity biomass at that timestep,
`w_t = X_{t\cdot\cdot} / \sum_t X_{t\cdot\cdot}`, and returns the weighted
sum reported as `BD_beta`.

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
