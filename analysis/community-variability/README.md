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
ifcb-process --dataset NESLTER_transect
ifcb-fill-missing data/NESLTER_transect
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

The analysis folder is workflow code rather than an installable package. To
fetch only the analysis and its metric-code dependencies:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set analysis/community-variability community.variability
```

Use clean data for observed-sample analyses:

```text
ifcb_carbon_clean.csv
ifcb_taxonomy.csv
```

Use fill data for balanced metacommunity analyses:

```text
ifcb_carbon_fill.csv
ifcb_taxonomy_fill.csv
```

## R Workflows

From this directory:

```bash
Rscript R/scripts/ifcb_single_season.R
Rscript R/scripts/ifcb_power_analysis.R
Rscript R/scripts/ifcb_sensitivity_analysis.R
Rscript R/scripts/ifcb_seasonal_comparison.R
```

## Python Workflows

Use fill data by default:

```bash
python python/scripts/ifcb_single_season.py --data-version fill
python python/scripts/ifcb_power_analysis.py --data-version fill
python python/scripts/ifcb_sensitivity_analysis.py --data-version fill
python python/scripts/ifcb_seasonal_comparison.py
```

Use clean data explicitly:

```bash
python python/scripts/ifcb_single_season.py --data-version clean
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
