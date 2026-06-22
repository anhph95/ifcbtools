# Python Analysis

The current directory is the analysis workspace, following the same movable
workspace convention used by `stingraytools`. By default, the scripts read:

```text
<current-directory>/data/NESLTER_transect
```

and write results and logs under:

```text
<current-directory>/results
```

Run the scripts from the intended workspace, or pass explicit `--data-dir`,
`--results-dir`, and `--log-dir` paths. Relative paths are resolved from the
current directory.

Install the Python metric and IFCB processing dependencies in the active
project environment. For example:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=community.variability/python"
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=ifcb.process"
```

Fetch only the Python analysis workflow and metric package:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set analysis/community-variability/python community.variability/python
```

```bash
python analysis/community-variability/python/scripts/ifcb_single_season.py --data-version fill
python analysis/community-variability/python/scripts/ifcb_power_analysis.py --data-version fill
python analysis/community-variability/python/scripts/ifcb_sensitivity_analysis.py --data-version fill
python analysis/community-variability/python/scripts/ifcb_seasonal_comparison.py
```

The Python workflow mirrors the R computational outputs and writes CSV products
for estimates, bootstrap summaries, sensitivity deltas, and seasonal comparison
tables.

Each script also writes timestamped `.out.log` and `.err.log` files under
`<results-dir>/logs`. Use `--log-level DEBUG` for more detail or
`--log-dir PATH` to choose another log directory; terminal logging remains
enabled. Each log begins with the command-line arguments and resolved workflow
paths; secret-like parameter names are automatically redacted.
