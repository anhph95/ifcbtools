# MATLAB Analysis

Open MATLAB in the folder where you want to run the analysis. The downloaded
scripts use that current folder as the analysis workspace: inputs are read from
`data/NESLTER_transect`, results are written under `results`, and logs are
written under `logs`.

Install or add the MATLAB metric functions separately, then download the
MATLAB analysis scripts into the current folder:

```matlab
eval(webread("https://raw.githubusercontent.com/anhph95/ifcbtools/main/analysis/community-variability/matlab/install_analysis_scripts.m"))
```

To replace existing copies, run:

```matlab
overwrite = true;
eval(webread("https://raw.githubusercontent.com/anhph95/ifcbtools/main/analysis/community-variability/matlab/install_analysis_scripts.m"))
```

Open the copied scripts and build on them. Keep `ifcb_common.m` in the same
folder as the analysis scripts because each script runs it with
`run("ifcb_common.m")`.

MATLAB command-window output and warnings are recorded with `diary` in
timestamped `.out.log` files under:

```text
logs/
```

Fatal exceptions are written with full reports to matching `.err.log` files.
Each output log begins with the data paths and analysis settings used for that
run.
