# MATLAB Analysis

Run these scripts from `analysis/community-variability` so `ifcb_common.m`
can resolve the repository paths and shared `results/` folder.

Add the MATLAB analysis and metric folders to the path from a checkout:

```matlab
run("analysis/community-variability/matlab/install_analysis_path.m")
install_analysis_path(true)
```

```matlab
run("matlab/scripts/ifcb_single_season.m")
run("matlab/scripts/ifcb_power_analysis.m")
run("matlab/scripts/ifcb_sensitivity_analysis.m")
run("matlab/scripts/ifcb_seasonal_comparison.m")
```

The MATLAB workflow reproduces the computational CSV outputs used by the R and
Python analysis scripts. The R workflow remains the plotting-rich reference for
publication-style figures.

MATLAB command-window output and warnings are recorded with `diary` in
timestamped `.out.log` files under:

```text
analysis/community-variability/logs/
```

Fatal exceptions are written with full reports to matching `.err.log` files.
Each output log begins with the data paths and analysis settings used for that
run.

To fetch only the MATLAB analysis and metric code from Git:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set analysis/community-variability/matlab community.variability/matlab
```
