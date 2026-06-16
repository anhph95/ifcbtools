# Community Variability

This folder contains the metacommunity variability code as separate MATLAB, R,
and Python implementations.

The shared ecological array convention is:

```text
X[time, site, taxon]
```

where `time` is repeated surveys, `site` is local communities, and `taxon` is
the IFCB taxon/class biomass dimension.

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
