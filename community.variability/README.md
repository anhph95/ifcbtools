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
metacommunity aggregate biomass. The notation below uses `\mu` for temporal
means, `\sigma` for temporal standard deviations, and `\sigma^2` for
variance terms.

### `CV_alpha`: Local Aggregate Variability

`CV_alpha` is temporal variability in aggregate biomass at the local-community
scale. It describes how much total biomass fluctuates within sites before
considering whether different sites fluctuate together.

$$
X_{ti\cdot} = \sum_j X_{tij},
\quad
\sigma_i = \mathrm{sd}_t(X_{ti\cdot}),
\quad
\mu_\gamma = \mathrm{mean}_t(X_{t\cdot\cdot})
$$

$$
CV_\alpha^2 =
\left(
\frac{\sum_i \sigma_i}
{\mu_\gamma}
\right)^2
$$

In code, summing the taxon dimension of `X[time, site, taxon]` gives
`X_{ti\cdot}` while retaining the site dimension. Each site therefore has its
own biomass time series. The implementation computes `\sigma_i` for each site,
sums those local standard deviations, and normalizes by `\mu_\gamma`, the mean
metacommunity biomass.

### `CV_gamma`: Metacommunity Aggregate Variability

`CV_gamma` is temporal variability in aggregate biomass at the metacommunity
scale. It describes how much total biomass fluctuates after all local
communities and taxa are pooled.

$$
X_{t\cdot\cdot} = \sum_i \sum_j X_{tij},
\quad
\mu_\gamma = \mathrm{mean}_t(X_{t\cdot\cdot}),
\quad
\sigma_\gamma = \mathrm{sd}_t(X_{t\cdot\cdot})
$$

$$
CV_\gamma^2 =
\left(
\frac{\sigma_\gamma}
{\mu_\gamma}
\right)^2
$$

In code, summing both the site and taxon dimensions produces one
metacommunity biomass time series, `X_{t\cdot\cdot}`. `CV_gamma()` computes
its mean `\mu_\gamma`, standard deviation `\sigma_\gamma`, and squared
coefficient of variation.

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

In code, `CV_phi` is computed as the ratio of the two stored aggregate
variability values. Because `CV_alpha` and `CV_gamma` are both squared
coefficients of variation, the ratio is the aggregate synchrony component of
the temporal alpha-gamma partition.

### `BD_alpha`: Local Compositional Variability

`BD_alpha` is temporal compositional variability at the local-community scale.
It describes how much taxon composition changes through time within sites,
after removing changes in total biomass.

$$
z_{tij} =
\sqrt{\frac{X_{tij}}{X_{ti\cdot}}},
\quad
\sigma^2_{ij} = \mathrm{Var}_t(z_{tij})
$$

$$
BD_i^h =
\sum_j \sigma^2_{ij}
$$

$$
w_i =
\frac{\mu_i}
{\sum_i \mu_i},
\quad
\mu_i = \mathrm{mean}_t(X_{ti\cdot})
$$

$$
BD_\alpha^h =
\sum_i w_i BD_i^h
$$

In code, local composition is represented by dividing each site-time taxon
vector by its local biomass and applying the Hellinger square-root
transformation. Temporal variance is then computed for each site-taxon
trajectory as `\sigma^2_{ij}`. Summing across taxa gives site-level
compositional variability, `BD_i^h`, and biomass weights `w_i` combine those
site-level values into `BD_alpha`.

### `BD_gamma`: Metacommunity Compositional Variability

`BD_gamma` is temporal compositional variability at the metacommunity scale.
It describes how regional taxon composition changes through time after pooling
local communities.

$$
X_{t\cdot j} = \sum_i X_{tij},
\quad
z_{t\cdot j} =
\sqrt{\frac{X_{t\cdot j}}{X_{t\cdot\cdot}}},
\quad
\sigma^2_{\gamma j} = \mathrm{Var}_t(z_{t\cdot j})
$$

$$
BD_\gamma^h =
\sum_j \sigma^2_{\gamma j}
$$

In code, each taxon is summed across sites to produce regional taxon biomass,
then divided by total metacommunity biomass to form regional composition.
`BD_gamma()` Hellinger-transforms that composition, computes temporal variance
for each regional taxon trajectory as `\sigma^2_{\gamma j}`, and sums those
variances across taxa.

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

In code, `BD_phi` is computed as the ratio of the two temporal compositional
variability values. The ratio is the compositional synchrony component of the
temporal alpha-gamma partition.

### `BD_beta`: Local-Scale Spatial Compositional Variability

`BD_beta` is compositional variability at the local spatial scale. It measures
how different local communities are from one another in composition at each
time, then summarizes that spatial variability through time.

$$
BD_t^h =
\sum_j \sigma^2_{tj},
\quad
\sigma^2_{tj} = \mathrm{Var}_i(z_{tij})
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
used here, but variance is computed across sites rather than through time. For
each timestep, `spatial_bd_by_time()` computes `\sigma^2_{tj}` for each taxon
and sums across taxa to obtain `BD_t^h`. `BD_spatial_weighted()` weights each
`BD_t^h` by total metacommunity biomass at that timestep and returns the
weighted sum reported as `BD_beta`.

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
