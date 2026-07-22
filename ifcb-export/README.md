# Python Processing Package

`ifcb` contains the NES-LTER IFCB processing workflow:

- MATLAB export from IFCB source products to local CSV files.
- Python cleaning, nearest-station assignment, bottle merge, nutrient merge.
- Optional fill product for analyses that need balanced station coverage.

## Install

From GitHub without manually cloning the repository:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=ifcb-export"
```

From a local checkout:

```bash
pip install -e ifcb-export
```

This installs:

```python
ifcb
```

and the commands:

```bash
ifcb
ifcb.process
```

To install the test dependency as well:

```bash
pip install -e "ifcb-export[test]"
```

## MATLAB Export

Set MATLAB's current folder to your workspace, then download the standalone
exporter:

```matlab
scriptUrl = "https://raw.githubusercontent.com/anhph95/ifcbtools/main/ifcb-export/matlab/export_ifcb.m";
websave("export_ifcb.m", scriptUrl);
```

Then run the downloaded function with one input MAT file:

```matlab
export_ifcb("\\sosiknas1\IFCB_products\NESLTER_transect\summary\carbon_group_class_withTS.mat")
```

With only the input path, the exporter scans available table variables, asks
you to confirm or choose the measurement table, infers whether the product is
`count` or `carbon` when possible, asks you to confirm before export, and asks
for an output folder relative to the current working directory. Press Enter to
use `data`, or enter a subfolder such as `data/NESLTER_transect`.
It processes one summary MAT file at a time and writes beneath the current
MATLAB workspace by default. The standard MAT filenames are:

```text
count_group_class_withTS.mat
carbon_group_class_withTS.mat
```

Scripted runs can provide the table and type explicitly. Carbon:

```matlab
export_ifcb("\\sosiknas1\IFCB_products\NESLTER_transect\summary\carbon_group_class_withTS.mat", data_table="classC_opt_adhoc_merge", data_type="carbon")
```

Count:

```matlab
export_ifcb("\\sosiknas1\IFCB_products\NESLTER_transect\summary\count_group_class_withTS.mat", data_table="classcount_opt_adhoc_merge", data_type="count")
```

Set `output_dir` explicitly for scripted runs. The function makes no
repository-location assumptions and does not search alternate paths.
The `data_table` and `data_type` arguments are optional for interactive MATLAB
use and recommended for scripted, reproducible runs.

## Data Products

The persisted processing products are:

```text
ifcb_count.csv / ifcb_carbon.csv  # MATLAB-exported measurement + metadata input
ifcb_class.csv                    # MATLAB-exported class/taxon lookup
*_clean*.csv                      # one selected Python output
*_fill.csv                        # optional analysis-specific fill output
```

## Processing Pipeline

`ifcb.process` always processes one explicitly selected input CSV.

Each pipeline operation is independently selectable and runs in the order
shown below, regardless of flag order:

```text
--clean -> --station -> --bottle -> --nutrient
```

Run the complete pipeline with:

```bash
ifcb data/NESLTER_transect/ifcb_count.csv --all
```

Unless `--output-file` is provided, `--all` writes:

```text
ifcb_count.csv -> ifcb_count_clean_station_bottle_nutrient.csv
```

`--clean` is the main workflow. It reads the selected MATLAB-exported
measurement + metadata CSV, runs quality-control checks, aggregates cast
replicates, and normalizes taxon values to per-liter floating-point values:

```bash
ifcb data/NESLTER_transect/ifcb_count.csv --clean --station --bottle --nutrient
```

When `--data-type` is omitted, the clean workflow infers `count` or `carbon`
from standard filenames such as `ifcb_count.csv` and `ifcb_carbon.csv`. Use
`--data-type` for ambiguous or renamed input files:

```bash
ifcb data/NESLTER_transect/ifcb_carbon.csv --clean --data-type carbon
```

To enrich an existing CSV without rerunning the clean workflow, omit
`--clean`:

```bash
ifcb data/NESLTER_transect/ifcb_count_clean.csv --station
ifcb data/NESLTER_transect/ifcb_count_clean.csv --bottle --nutrient
```

At least one processing operation must be selected.

Unless `--output-file` is given, selected operation names are appended before
the input extension:

```text
sample.csv --clean                         -> sample_clean.csv
sample.csv --station                   -> sample_station.csv
sample.csv --clean --station           -> sample_clean_station.csv
sample.csv --bottle --nutrient -> sample_bottle_nutrient.csv
sample.csv --all                           -> sample_clean_station_bottle_nutrient.csv
```

The pipeline writes one final CSV for the selected operations. Intermediate
dataframes are kept in memory.

During cleaning, `ifcb_taxonomy.csv` defaults to a file beside the input CSV.
Override that path when needed:

```bash
ifcb data/NESLTER_transect/counts.csv --clean --data-type count --taxonomy-file taxa.csv --output-file counts_clean.csv
```

The `--clean` operation runs these steps in order:

1. Read one MATLAB-exported count or carbon CSV that already contains metadata.
2. Require metadata columns `skip`, `sample_type`, `depth`, `ml_analyzed`,
   `longitude`, `latitude`, and `sample_time`.
3. Keep only rows where `skip == 0`.
4. Convert `depth`, `ml_analyzed`, `longitude`, and `latitude` to numeric
   values.
5. For `underway`, `underway_discrete`, and `bucket` rows, fill missing
   `depth` with `0`.
6. Replace `ml_analyzed == 0` with missing.
7. Drop rows missing `sample_time`, `longitude`, `latitude`, or
   `ml_analyzed`.
8. Parse `sample_time` as a timestamp and drop rows with unparseable times.
9. Add `year`, `month`, `day`, ISO `week`, day-of-year `doy`, and seasonal
    quarter `season`.
10. Add empty `nearest_station` and `station_distance` columns for later
    station assignment.
11. Read `ifcb_taxonomy.csv` and use its `Annotations` column to select taxon
    columns from the count or carbon input table.
12. Aggregate `cast` replicate rows by `cruise`, `cast`, and `depth`.
    Taxon columns, `ml_analyzed`, and `n_images` are summed; `sample_time` is
    the earliest replicate time; other metadata fields keep one unique value
    or join multiple unique values.
13. Resolve the data type from `--data-type` or from standard input filenames
    containing `count` or `carbon`.
14. Normalize each selected taxon column to product-specific per-liter units:
    count uses `normalized_value = raw_value / ml_analyzed * 1000` for
    cells L-1, while carbon uses
    `normalized_value = raw_value / ml_analyzed * 0.001` for ug C L-1.

After `--clean`, the optional pipeline steps run in this order when selected:

1. `--station` assigns nearest-station fields.
2. `--bottle` merges CTD bottle fields.
3. `--nutrient` merges nutrient fields.

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
python -m ifcb.fill --input data/NESLTER_transect/ifcb_carbon_clean_station_bottle_nutrient.csv
```

This reads one explicit input CSV that already includes `nearest_station` from
the station-assignment step, creates missing `cast_from_udw` rows from
same-cruise underway samples with existing station assignments, fills
bottle/nutrient values only for the new rows, and writes by default:

```text
ifcb_carbon_clean_station_bottle_nutrient.csv -> ifcb_carbon_clean_station_bottle_nutrient_fill.csv
```

Rows that still lack `nearest_station` are ignored as fill candidates; run the
station step before fill.

Input and output paths can also be selected explicitly:

```bash
python -m ifcb.fill --input data/NESLTER_transect/biomass_cleaned.csv --output-file biomass_filled.csv
```

## Taxonomy Mapping

Use `ifcb.taxonomy` when taxon columns need to be aggregated from one taxonomy
level to another.

```bash
ifcb.taxonomy --input data/NESLTER_transect/ifcb_carbon_mix.csv --taxonomy-file data/NESLTER_transect/ifcb_taxonomy.csv --output data/NESLTER_transect/ifcb_carbon_label.csv --from-level Annotations --to-level Label
```

The same mapper is available in Python:

```python
from ifcb.taxonomy import taxonomy_mapping
```

## Useful Imports

```python
from ifcb import process
from ifcb.fill import fill_dataset
from ifcb.stations import StationLocator, nearest_station
```
