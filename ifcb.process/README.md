# Python Processing Package

`ifcb.process` contains the NES-LTER IFCB processing workflow:

- MATLAB export from IFCB source products to local CSV files.
- Python cleaning, nearest-station assignment, bottle merge, nutrient merge.
- Optional fill product for analyses that need balanced station coverage.

## Install

From GitHub without manually cloning the repository:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=ifcb.process"
```

From a local checkout:

```bash
pip install -e ifcb.process
```

This installs:

```python
ifcb.process.neslter
```

and the commands:

```bash
ifcb-process
ifcb-fill-missing
```

## MATLAB Export

From MATLAB, run:

```matlab
run("ifcb.process/matlab/export_ifcb_mat.m")
```

To add the export folder to the MATLAB path:

```matlab
run("ifcb.process/matlab/install_ifcb_process_export.m")
install_ifcb_process_export(true)
```

To fetch only this package from Git:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set ifcb.process
```

## Data Products

The persisted processing products are:

```text
*_raw.csv    # MATLAB-exported input
*_clean.csv  # cleaned, station-assigned, bottle-merged, nutrient-merged
*_fill.csv   # optional fill product for balanced metacommunity analyses
```

## Clean Processing

```bash
ifcb-process --dataset NESLTER_broadscale
ifcb-process data/NESLTER_transect
```

Expected inputs in the dataset directory:

```text
ifcb_metadata.csv
ifcb_taxonomy.csv
ifcb_count_raw.csv
ifcb_carbon_raw.csv
```

`ifcb-process` writes:

```text
ifcb_count_clean.csv
ifcb_carbon_clean.csv
```

Clean data include:

- metadata cleaning
- cast replicate aggregation
- count/carbon normalization
- nearest station assignment
- CTD bottle merge
- nutrient merge

## Optional Fill Product

Use this only when an analysis needs balanced station coverage.

```bash
ifcb-fill-missing data/NESLTER_transect
```

This reads `*_clean.csv`, creates missing `cast_from_udw` rows from same-cruise
underway samples, fills bottle/nutrient values only for the new rows, and writes:

```text
ifcb_count_fill.csv
ifcb_carbon_fill.csv
ifcb_taxonomy_fill.csv
```

## Useful Imports

```python
from ifcb.process.neslter.process import process
from ifcb.process.neslter.fill import make_filled_dataset
from ifcb.process.neslter.stations import StationLocator, nearest_station
```
