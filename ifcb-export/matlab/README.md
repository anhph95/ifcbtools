# MATLAB Export

This folder exports one source IFCB MATLAB product to CSV files used by the
Python processing package.

## Download and Run

Set MATLAB's current folder to the workspace where the `data/` directory
should be created. Download the standalone script directly from GitHub:

```matlab
scriptUrl = "https://raw.githubusercontent.com/anhph95/ifcbtools/main/ifcb-export/matlab/export_ifcb.m";
websave("export_ifcb.m", scriptUrl);
```

Then run the downloaded function with one input MAT file:

```matlab
export_ifcb( ...
    "\\sosiknas1\IFCB_products\NESLTER_transect\summary\carbon_group_class_withTS.mat")
```

With only the input path, the function scans available table variables, asks
you to confirm or choose the measurement table, infers whether the product is
`count` or `carbon` when possible, asks you to confirm before export, and asks
for an output folder relative to the current working directory. Press Enter to
use `data`, or enter a subfolder such as `data/NESLTER_transect`.
The function does not depend on the repository layout. It processes one MAT
file at a time. Standard inputs are:

```text
count_group_class_withTS.mat
carbon_group_class_withTS.mat
```

For scripted carbon export, provide the table and type explicitly:

```matlab
export_ifcb( ...
    "\\sosiknas1\IFCB_products\NESLTER_transect\summary\carbon_group_class_withTS.mat", ...
    data_table="classC_opt_adhoc_merge", ...
    data_type="carbon")
```

For scripted count export, select the count MAT file and count data table:

```matlab
export_ifcb( ...
    "\\sosiknas1\IFCB_products\NESLTER_transect\summary\count_group_class_withTS.mat", ...
    data_table="classcount_opt_adhoc_merge", ...
    data_type="count")
```

For scripted runs, set the output path explicitly when needed:

```matlab
export_ifcb( ...
    "path/to/product.mat", ...
    data_table="custom_carbon_table", ...
    data_type="carbon", ...
    output_dir="path/to/output")
```

The function does not search for alternate paths.

## Expected Output

When prompted for an output folder, press Enter to write local files under:

```text
<workdir>/data/
```

Command-window output is also written to a timestamped `.out.log`, and fatal
exceptions to a matching `.err.log`, under:

```text
<workdir>/logs/
```

This is the same default logging folder used by the Python processing commands
when they are invoked from the same working directory.
The output log records the selected dataset, source directory, input files, and
resolved output and logging directories.

The Python processing step expects:

```text
ifcb_count.csv or ifcb_carbon.csv
ifcb_class.csv
ifcb_taxonomy.csv
```

The count or carbon CSV contains measurement columns merged with `meta_data` by
`pid`. `ifcb_class.csv` is exported separately as the class/taxon lookup and is
not merged into the measurement table. If `ifcb_taxonomy.csv` is missing,
`ifcb.process` can download the configured Google Sheet taxonomy during Python
processing.
