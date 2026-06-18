# Python Analysis

Run these scripts from `analysis/community-variability` or pass explicit
`--data-dir` and `--results-dir` paths.

Install the Python metric dependency directly from Git:

```bash
pip install "git+https://github.com/anhph95/ifcbtools.git#subdirectory=community.variability/python"
```

Fetch only the Python analysis workflow and metric package:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set analysis/community-variability/python community.variability/python
```

```bash
python python/scripts/ifcb_single_season.py --data-version fill
python python/scripts/ifcb_power_analysis.py --data-version fill
python python/scripts/ifcb_sensitivity_analysis.py --data-version fill
python python/scripts/ifcb_seasonal_comparison.py
```

The Python workflow mirrors the R computational outputs and writes CSV products
for estimates, bootstrap summaries, sensitivity deltas, and seasonal comparison
tables.

Each script also writes timestamped `.out.log` and `.err.log` files under
`<results-dir>/logs`. Use `--log-level DEBUG` for more detail or
`--log-dir PATH` to choose another log directory; terminal logging remains
enabled. Each log begins with the command-line arguments and resolved workflow
paths; secret-like parameter names are automatically redacted.
