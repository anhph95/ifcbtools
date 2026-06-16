# ifcbtools

Mixed Python and R tools for NES-LTER IFCB data processing and community
variability analysis.

## Repository layout

```text
ifcbtools/
|-- pyproject.toml
|-- scripts/matlab/                   # MATLAB export step
|-- data/                             # Local MATLAB exports, ignored by Git
|-- src/ifcb/neslter/                 # Python data processing package
|-- R/community.variability/          # R package
`-- analysis/community-variability/   # R workflows and generated results
```

## Install the Python package

```bash
pip install -e .
```

This installs the `ifcb.neslter` package and the `ifcb-process` command.

## Workflow

```text
MATLAB export -> Python data process -> R analysis
```

## Export MATLAB products

Run the MATLAB export script:

```matlab
run("scripts/matlab/export_ifcb_mat.m")
```

This writes local CSV files to:

```text
data/<dataset>/
```

## Install the R package

```bash
R CMD INSTALL R/community.variability
```

or from an R session:

```r
install.packages("R/community.variability", repos = NULL, type = "source")
```

The R package provides metacommunity variability functions for arrays:

```r
X[t, i, j]
```

where `t` is timestep, `i` is station or local community, and `j` is taxon.

## Run the Python CLI

```bash
ifcb-process
ifcb-process --dataset NESLTER_broadscale
ifcb-process /path/to/input_data -o /path/to/output --sample-type cast underway
```

Expected MATLAB-exported input files:

- `ifcb_metadata.csv`
- `ifcb_class.csv`
- `ifcb_count_raw.csv`
- `ifcb_carbon_raw.csv`

The Python process uses `ifcb_taxonomy.csv` to identify taxon columns. If that
file is missing, `ifcb-process` downloads it from the configured Google Sheet
and saves it beside the MATLAB exports.

## Use station lookup independently

```python
from ifcb.neslter.stations import nearest_station, assign_nearest_stations

station, distance_km = nearest_station(
    lat=41.2,
    lon=-70.8,
    timestamp="2023-07-01T12:00:00",
    station_reference="stations.csv",  # optional; defaults to NES-LTER API
)

with_stations = assign_nearest_stations(samples_df, station_reference="stations.csv")
```

## Run the R workflows

Install the R package first, then run scripts from the analysis project
directory:

```bash
cd analysis/community-variability
Rscript scripts/ifcb_single_season.R
Rscript scripts/ifcb_power_analysis.R
Rscript scripts/ifcb_sensitivity_analysis.R
Rscript scripts/ifcb_seasonal_comparison.R
```
