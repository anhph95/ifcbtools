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

Set MATLAB's current folder to your workspace, then download and run the
standalone exporter:

```matlab
scriptUrl = "https://raw.githubusercontent.com/anhph95/ifcbtools/main/ifcb.process/matlab/export_ifcb.m";
websave("export_ifcb.m", scriptUrl);
export_ifcb()
```

The default call processes `NESLTER_transect` from the server and writes
beneath the current MATLAB workspace. The standard MAT filenames are:

```text
count_group_class_withTS.mat
carbon_group_class_withTS.mat
```

Run another dataset without editing the function:

```matlab
export_ifcb("NESLTER_broadscale")
```

Override `SummaryDir`, `CountMatFile`, `CarbonMatFile`, or `OutputDir` only
when needed. The function makes no repository-location assumptions and does
not search alternate paths.

## Data Products

The persisted processing products are:

```text
*_raw.csv    # MATLAB-exported input
*_clean.csv  # selected output; contents depend on the requested operations
*_fill.csv   # optional fill product for balanced metacommunity analyses
```

## Processing Pipeline

`ifcb-process` always processes one explicitly selected input CSV.

Each pipeline operation is independently selectable and runs in the order
shown below, regardless of flag order:

```text
--clean -> --add-station -> --merge-bottle -> --merge-nutrient
```

`--clean` is the main workflow. It reads one raw count or carbon CSV, filters
and cleans metadata, aggregates cast replicates, and normalizes taxon values.
Select the input type explicitly:

```bash
ifcb-process data/NESLTER_transect/ifcb_count_raw.csv \
  --clean \
  --data-type count \
  --add-station \
  --merge-bottle \
  --merge-nutrient
```

To enrich an existing CSV without rerunning the clean workflow, omit
`--clean`:

```bash
ifcb-process data/NESLTER_transect/ifcb_count_clean.csv --add-station
ifcb-process data/NESLTER_transect/ifcb_count_clean.csv --merge-bottle --merge-nutrient
```

At least one processing operation must be selected.

Unless `--output-file` is given, selected operation names are appended before
the input extension:

```text
sample.csv --clean                         -> sample_clean.csv
sample.csv --add-station                   -> sample_station.csv
sample.csv --clean --add-station           -> sample_clean_station.csv
sample.csv --merge-bottle --merge-nutrient -> sample_bottle_nutrient.csv
```

During cleaning, `ifcb_metadata.csv` and `ifcb_taxonomy.csv` default to files
beside the input CSV. Override either path when needed:

```bash
ifcb-process data/NESLTER_transect/counts.csv \
  --clean \
  --data-type count \
  --metadata-file metadata.csv \
  --taxonomy-file taxa.csv \
  --output-file counts_processed.csv
```

The complete pipeline includes:

- metadata cleaning
- cast replicate aggregation
- count/carbon normalization
- nearest station assignment
- CTD bottle merge
- nutrient merge

## Logging

Each Python processing command records the workflow steps in timestamped files:

```text
<output-dir>/logs/<command>_YYYYMMDD_HHMMSS.out.log
<output-dir>/logs/<command>_YYYYMMDD_HHMMSS.err.log
```

The `.out.log` file records messages at the selected level and above, while
the `.err.log` file records errors. The same messages remain visible in the
terminal. Use `--log-level DEBUG` for more detail or `--log-dir PATH` to place
the files elsewhere. Each log begins with the resolved command-line inputs and
paths used for that run; secret-like parameter names are automatically redacted.

## Optional Fill Product

Use this only when an analysis needs balanced station coverage.

```bash
ifcb-fill-missing data/NESLTER_transect
```

This reads the default `*_clean.csv` files, creates missing `cast_from_udw` rows from same-cruise
underway samples, fills bottle/nutrient values only for the new rows, and writes:

```text
ifcb_count_fill.csv
ifcb_carbon_fill.csv
```

The taxonomy CSV is read to map annotations into `Label` groups but is not
modified or copied. Input and output filenames can also be selected explicitly:

```bash
ifcb-fill-missing data/NESLTER_transect \
  --taxonomy-file taxa.csv \
  --count-file counts_cleaned.csv \
  --carbon-file biomass_cleaned.csv \
  --count-output-file counts_filled.csv \
  --carbon-output-file biomass_filled.csv
```

## Useful Imports

```python
from ifcb.process.neslter.process import process
from ifcb.process.neslter.fill import make_filled_dataset
from ifcb.process.neslter.stations import StationLocator, nearest_station
```
