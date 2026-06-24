# ifcbtools

Mixed MATLAB, Python, R, and analysis tools for NES-LTER IFCB data processing
and metacommunity variability analysis.

## Repository Layout

```text
ifcbtools/
|-- ifcb.process/                   # MATLAB export plus Python processing package
|-- community.variability/          # Installable MATLAB, R, and Python metric code
|-- analysis/                       # Analysis workflows using the packages above
|-- data/                           # Local data products, ignored by Git
`-- README.md
```

## Main Workflow

```text
ifcb.process MATLAB export -> ifcb.process Python clean -> optional fill -> community.variability analysis
```

## Documentation

- [IFCB processing](ifcb.process/README.md)
- [Community variability](community.variability/README.md)
- [Community variability analysis](analysis/community-variability/README.md)

## Required MATLAB Export

All downstream processing starts from CSV files exported from the IFCB MATLAB
summary products. Set MATLAB's current folder to your workspace, then download
and run the standalone exporter:

```matlab
scriptUrl = "https://raw.githubusercontent.com/anhph95/ifcbtools/main/ifcb.process/matlab/export_ifcb.m";
websave("export_ifcb.m", scriptUrl);
export_ifcb()
```

By default, the exporter processes `NESLTER_transect`. It reads from:

```text
\\sosiknas1\IFCB_products\<dataset>\summary
```

and write the required exported products under the active MATLAB workspace:

```text
<workdir>/data/<dataset>/ifcb_class.csv
<workdir>/data/<dataset>/ifcb_metadata.csv
<workdir>/data/<dataset>/ifcb_count.csv
<workdir>/data/<dataset>/ifcb_carbon.csv
```

Process a different dataset without editing the downloaded function:

```matlab
export_ifcb("NESLTER_broadscale")
```

The standard count and carbon MAT filenames remain the defaults. Override
`SummaryDir`, `CountMatFile`, `CarbonMatFile`, or `OutputDir` only when needed.
The function does not assume a repository layout or search alternate paths.

## Install Processing Package

From GitHub without manually cloning the repository:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=ifcb.process"
```

From a local checkout:

```bash
pip install -e ifcb.process
```

This installs the `ifcb.process` package plus the processing commands:

```bash
ifcb-process
ifcb-fill-missing
```

## Install Community Variability Packages

Python from GitHub:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=community.variability/python"
```

R from GitHub:

```r
install.packages("remotes")
remotes::install_git(
  "https://github.com/anhph95/ifcbtools.git",
  subdir = "community.variability/R/community.variability"
)
```

Local checkout:

```bash
pip install -e community.variability/python
R CMD INSTALL community.variability/R/community.variability
```

For MATLAB, add the metric folder to the MATLAB path:

```matlab
addpath("community.variability/matlab/community_variability")
```

or run from a checkout:

```matlab
run("community.variability/matlab/install_community_variability.m")
install_community_variability(true)
```

## Fetch Analysis Only

Install the metric package separately, then use Git sparse checkout only to
copy the analysis workflow:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set analysis/community-variability
```
