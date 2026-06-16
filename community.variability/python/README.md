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
