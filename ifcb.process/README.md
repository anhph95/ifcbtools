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
ifcb_count.csv / ifcb_carbon.csv  # MATLAB-exported input
*_clean.csv                       # selected output; contents depend on the requested operations
*_fill.csv                        # optional fill product for balanced metacommunity analyses
```

## Processing Pipeline

`ifcb-process` always processes one explicitly selected input CSV.

Each pipeline operation is independently selectable and runs in the order
shown below, regardless of flag order:

```text
--clean -> --add-station -> --merge-bottle -> --merge-nutrient
```

Run the complete pipeline with:

```bash
ifcb-process data/NESLTER_transect/ifcb_count.csv --all
```

Unless `--output-file` is provided, `--all` writes:

```text
ifcb_count.csv -> ifcb_count_clean_station_bottle_nutrient.csv
```

`--clean` is the main workflow. It reads the selected MATLAB-exported CSV,
cleans metadata, aggregates cast replicates, and normalizes taxon values to
per-liter floating-point values:

```bash
ifcb-process data/NESLTER_transect/ifcb_count.csv \
  --clean \
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
sample.csv --all                           -> sample_clean_station_bottle_nutrient.csv
```

During cleaning, `ifcb_metadata.csv` and `ifcb_taxonomy.csv` default to files
beside the input CSV. Override either path when needed:

```bash
ifcb-process data/NESLTER_transect/counts.csv \
  --clean \
  --metadata-file metadata.csv \
  --taxonomy-file taxa.csv \
  --output-file counts_clean.csv
```

The `--clean` operation runs these steps in order:

1. Read `ifcb_metadata.csv`.
2. Require metadata columns `skip`, `sample_type`, `depth`, `ml_analyzed`,
   `longitude`, `latitude`, and `sample_time`.
3. Keep only rows where `skip == 0`.
4. Keep only `sample_type` values `cast`, `underway`, and
   `underway_discrete`. Override this set with `--sample-type`.
5. Rewrite `sample_type == "underway_discrete"` to `sample_type == "underway"`.
6. Convert `depth`, `ml_analyzed`, `longitude`, and `latitude` to numeric
   values.
7. For `underway` rows, fill missing `depth` with `0`.
8. Replace `ml_analyzed == 0` with missing.
9. Drop rows missing `sample_time`, `longitude`, `latitude`, or
   `ml_analyzed`.
10. Parse `sample_time` as a timestamp and drop rows with unparseable times.
11. Add `year`, `month`, `day`, ISO `week`, day-of-year `doy`, and seasonal
    quarter `season`.
12. Add empty `nearest_station` and `station_distance` columns for later
    station assignment.
13. Read `ifcb_taxonomy.csv` and use its `Annotations` column to select taxon
    columns from the count or carbon input table.
14. Merge cleaned metadata and the selected count or carbon table by `pid`.
15. Aggregate `cast` replicate rows by `cruise`, `cast`, and `depth`.
    Taxon columns, `ml_analyzed`, and `n_images` are summed; `sample_time` is
    the earliest replicate time; other metadata fields keep one unique value
    or join multiple unique values.
16. Normalize each selected taxon column to product-specific per-liter units:
    count files use `normalized_value = raw_value / ml_analyzed * 1000`
    for cells L-1, while carbon files use
    `normalized_value = raw_value / ml_analyzed * 0.001` for ug C L-1.

After `--clean`, the optional pipeline steps run in this order when selected:

1. `--add-station` assigns nearest-station fields.
2. `--merge-bottle` merges CTD bottle fields.
3. `--merge-nutrient` merges nutrient fields.

## Logging

Each Python processing command records the workflow steps in timestamped files:

```text
<current-directory>/logs/<command>_YYYYMMDD_HHMMSS.out.log
<current-directory>/logs/<command>_YYYYMMDD_HHMMSS.err.log
```

The `.out.log` file records messages at the selected level and above, while
the `.err.log` file records errors. The same messages remain visible in the
terminal. Use `--log-level DEBUG` for more detail or `--log-dir PATH` to place
the files elsewhere. Each log begins with the resolved command-line inputs and
paths used for that run; secret-like parameter names are automatically redacted.

## Optional Fill Product

Use this only when an analysis needs balanced station coverage.

```bash
ifcb-fill-missing data/NESLTER_transect/ifcb_carbon_clean.csv
```

This reads one explicit input CSV, creates missing `cast_from_udw` rows from
same-cruise underway samples, fills bottle/nutrient values only for the new
rows, and writes by default:

```text
ifcb_carbon_clean.csv -> ifcb_carbon_clean_fill.csv
```

The taxonomy CSV is read to map annotations into `Label` groups but is not
modified or copied. If `nearest_station` is absent, nearest stations are
assigned automatically before filling.

Input and output paths can also be selected explicitly:

```bash
ifcb-fill-missing data/NESLTER_transect/biomass_cleaned.csv \
  --taxonomy-file taxa.csv \
  --output-file biomass_filled.csv
```

## Useful Imports

```python
from ifcb.process.neslter.process import process
from ifcb.process.neslter.fill import make_filled_dataset
from ifcb.process.neslter.stations import StationLocator, nearest_station
```
