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

## Use

```python
from community_variability import make_community_array, calc_metacommunity_metrics

community = make_community_array(data_wide, taxa_cols)
metrics = calc_metacommunity_metrics(community)
```

The core array convention is:

```text
X[time, site, taxon]
```

## Mathematical Notation

The Python implementation uses the same notation as the top-level community
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
`make_community_array()` builds `X[time, site, taxon]`. Metric functions use
NumPy dimension reductions such as `np.nansum()`, `np.nanstd()`, and
`np.nanvar()` so the code follows the equations directly.

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
