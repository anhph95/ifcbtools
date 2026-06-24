# R Analysis

Run these scripts from `analysis/community-variability` so all outputs land in
the shared `results/` folder.

Install the R metric dependency directly from Git:

```r
install.packages("remotes")
remotes::install_git(
  "https://github.com/anhph95/ifcbtools.git",
  subdir = "community.variability/R/community.variability"
)
```

Fetch only the R analysis workflow and R metric package:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set analysis/community-variability/R community.variability/R
```

```bash
Rscript R/scripts/ifcb_single_season.R
Rscript R/scripts/ifcb_power_analysis.R
Rscript R/scripts/ifcb_sensitivity_analysis.R
Rscript R/scripts/ifcb_seasonal_comparison.R
```

The R scripts are the plotting-rich reference workflow for the analysis.

Each workflow writes timestamped `.out.log` and `.err.log` files under:

```text
analysis/community-variability/logs/
```

Informational messages remain visible in the terminal. Fatal errors are also
written to the error log with traceback information. Each output log records
the workflow settings and paths used for that run.
