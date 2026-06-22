# MATLAB Export

This folder exports source IFCB MATLAB products to CSV files used by the Python
processing package.

## Download and Run

Set MATLAB's current folder to the workspace where the `data/` directory
should be created. Download the standalone script directly from GitHub:

```matlab
scriptUrl = "https://raw.githubusercontent.com/anhph95/ifcbtools/main/ifcb.process/matlab/export_ifcb.m";
websave("export_ifcb.m", scriptUrl);
run("export_ifcb.m")
```

The script does not depend on the repository layout. Its editable user
settings keep the NES-LTER server as the default source:

```matlab
dataset = 'NESLTER_broadscale';
summaryDir = fullfile('\\sosiknas1', 'IFCB_products', dataset, 'summary');
outputDir = fullfile(pwd, 'data', dataset);
```

Change `summaryDir` or `outputDir` directly when the source MAT products or
destination are elsewhere. The script does not search for alternate paths.

## Expected Output

By default, the export writes local files under the MATLAB working directory:

```text
<workdir>/data/<dataset>/
```

Command-window output is also written to a timestamped `.out.log`, and fatal
exceptions to a matching `.err.log`, under:

```text
data/<dataset>/logs/
```

This is the same logging folder used by the Python processing commands.
The output log records the selected dataset, source directory, input files, and
resolved output and logging directories.

The Python processing step expects:

```text
ifcb_metadata.csv
ifcb_count_raw.csv
ifcb_carbon_raw.csv
ifcb_taxonomy.csv
```

If `ifcb_taxonomy.csv` is missing, `ifcb-process` can download the configured
Google Sheet taxonomy during Python processing.
