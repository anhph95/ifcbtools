# MATLAB Export

This folder exports source IFCB MATLAB products to CSV files used by the Python
processing package.

## Run Export

From MATLAB, run:

```matlab
run("ifcb.process/matlab/export_ifcb_mat.m")
```

To add this folder to the MATLAB path:

```matlab
run("ifcb.process/matlab/install_ifcb_process_export.m")
install_ifcb_process_export(true)
```

## Fetch Only This Folder From Git

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set ifcb.process/matlab
```

## Expected Output

The export writes local files under:

```text
data/<dataset>/
```

The Python processing step expects:

```text
ifcb_metadata.csv
ifcb_count_raw.csv
ifcb_carbon_raw.csv
ifcb_taxonomy.csv
```

If `ifcb_taxonomy.csv` is missing, `ifcb-process` can download the configured
Google Sheet taxonomy during Python processing.
