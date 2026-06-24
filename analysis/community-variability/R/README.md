# R Analysis

Open R or RStudio in the directory where you want to run the analysis. The
downloaded scripts use that current working directory as the analysis workspace:
inputs are read from `data/NESLTER_transect`, results are written under
`results`, and logs are written under `logs`.

Install the R metric dependency directly from Git:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_git(
  "https://github.com/anhph95/ifcbtools.git",
  subdir = "community.variability/R/community.variability"
)
```

Download the R analysis scripts into the current working directory:

```r
source("https://raw.githubusercontent.com/anhph95/ifcbtools/main/analysis/community-variability/R/install_analysis_scripts.R")
```

To replace existing copies, run:

```r
options(ifcb.analysis.overwrite = TRUE)
source("https://raw.githubusercontent.com/anhph95/ifcbtools/main/analysis/community-variability/R/install_analysis_scripts.R")
```

Open the copied scripts and run them interactively. Keep `ifcb_common.R` in the
same directory as the analysis scripts because each script loads it with
`source("ifcb_common.R")`.

Each workflow writes timestamped `.out.log` and `.err.log` files under:

```text
logs/
```

Informational messages remain visible in the terminal. Fatal errors are also
written to the error log with traceback information. Each output log records
the workflow settings and paths used for that run.
