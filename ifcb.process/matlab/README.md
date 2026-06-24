# MATLAB Export

This folder exports source IFCB MATLAB products to CSV files used by the Python
processing package.

## Download and Run

Set MATLAB's current folder to the workspace where the `data/` directory
should be created. Download the standalone script directly from GitHub:

```matlab
scriptUrl = "https://raw.githubusercontent.com/anhph95/ifcbtools/main/ifcb.process/matlab/export_ifcb.m";
websave("export_ifcb.m", scriptUrl);
export_ifcb()
```

The function does not depend on the repository layout. By default it processes
`NESLTER_transect` from the server using:

```text
count_group_class_withTS.mat
carbon_group_class_withTS.mat
```

Run another dataset without editing the function:

```matlab
export_ifcb("NESLTER_broadscale")
```

Override paths only when needed:

```matlab
export_ifcb("NESLTER_broadscale", ...
    SummaryDir="path/to/summary", ...
    OutputDir="path/to/output")
```

`CountMatFile` and `CarbonMatFile` can also override either standard filename.
The function does not search for alternate paths.

## Expected Output

By default, the export writes local files under the MATLAB working directory:

```text
<workdir>/data/<dataset>/
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
ifcb_metadata.csv
ifcb_count.csv
ifcb_carbon.csv
ifcb_taxonomy.csv
```

If `ifcb_taxonomy.csv` is missing, `ifcb-process` can download the configured
Google Sheet taxonomy during Python processing.
