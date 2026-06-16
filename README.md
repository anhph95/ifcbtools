# ifcbtools

Mixed Python and R tools for NES-LTER IFCB data processing and community
variability analysis.

## Repository layout

```text
ifcbtools/
|-- pyproject.toml
|-- data/                             # Local MATLAB exports, ignored by Git
|-- matlab/                           # MATLAB export step
|-- python/ifcb/neslter/              # Python data processing package
|-- r/community.variability/          # R package
`-- analysis/community-variability/   # R workflows and generated results
```

## Install the Python package

```bash
pip install -e .
```

This installs the `ifcb.neslter` package and the staged processing commands.

## Workflow

```text
MATLAB export
  -> clean
  -> add nearest_station
  -> optionally fill missing cast data
  -> merge bottle data
  -> merge nutrient data
  -> analysis
```

## Export MATLAB products

Run the MATLAB export script:

```matlab
run("matlab/export_ifcb_mat.m")
```

This writes local CSV files to:

```text
data/<dataset>/
```

## Install the R package

```bash
R CMD INSTALL r/community.variability
```

or from an R session:

```r
install.packages("r/community.variability", repos = NULL, type = "source")
```

The R package provides metacommunity variability functions for arrays:

```r
X[t, i, j]
```

where `t` is timestep, `i` is station or local community, and `j` is taxon.

## Run the Python pipeline

The Python pipeline is split into explicit stages so different analyses can use
different products.

### 1. Clean MATLAB exports

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

This stage writes:

```text
ifcb_count_clean.csv
ifcb_carbon_clean.csv
```

### 2. Add nearest station

Assign `nearest_station` and `station_distance` to all rows using the current
`StationLocator` workflow:

```bash
ifcb-add-stations data/NESLTER_transect
```

This reads `*_clean.csv` and writes:

```text
ifcb_count_station.csv
ifcb_carbon_station.csv
```

### 3. Optionally fill missing cast data

For balanced metacommunity analyses, fill missing surface or station casts from
same-cruise underway samples:

```bash
ifcb-fill-missing data/NESLTER_transect --data-type carbon
```

This reads `*_station.csv` and writes:

```text
ifcb_carbon_filled.csv
ifcb_taxonomy_filled.csv
```

The filled product is optional. Use `*_clean.csv` or `*_station.csv` for
analyses that should only use observed samples.

### 4. Merge bottle data

```bash
ifcb-merge-bottle data/NESLTER_transect --data-type carbon
```

By default this reads `*_filled.csv` and writes:

```text
ifcb_carbon_bottle.csv
```

Use `--input-stage station` or another stage name if you want bottle data
merged into a different product.

### 5. Merge nutrient data

```bash
ifcb-merge-nutrient data/NESLTER_transect --data-type carbon
```

By default this reads `*_bottle.csv` and writes:

```text
ifcb_carbon_nutrient.csv
```

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

## Run the Python community workflows

The Python community-variability scripts can use either clean or filled data.
For balanced seasonal metacommunity analyses, use the filled product:

```bash
python analysis/community-variability/scripts/ifcb_single_season.py --data-version filled
python analysis/community-variability/scripts/ifcb_power_analysis.py --data-version filled
python analysis/community-variability/scripts/ifcb_sensitivity_analysis.py --data-version filled
python analysis/community-variability/scripts/ifcb_seasonal_comparison.py
```

To use only the clean product:

```bash
python analysis/community-variability/scripts/ifcb_single_season.py --data-version clean
```
