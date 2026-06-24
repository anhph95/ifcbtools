# Python Analysis

Open Python in the directory where you want to run the analysis. The downloaded
scripts use that current working directory as the analysis workspace. By
default, the scripts read:

```text
<current-directory>/data/NESLTER_transect
```

write results under:

```text
<current-directory>/results
```

By default, logs are written under:

```text
<current-directory>/logs
```

Run the scripts from the intended workspace, or pass explicit `--data-dir`,
`--results-dir`, and `--log-dir` paths. Relative paths are resolved from the
current directory.

Install the Python metric and IFCB processing dependencies separately in the
active project environment.

Download the Python analysis scripts into the current working directory:

```python
import urllib.request

exec(
    urllib.request.urlopen(
        "https://raw.githubusercontent.com/anhph95/ifcbtools/main/analysis/community-variability/python/install_analysis_scripts.py"
    )
    .read()
    .decode("utf-8")
)
```

To replace existing copies, run:

```python
import urllib.request

IFCB_ANALYSIS_OVERWRITE = True
exec(
    urllib.request.urlopen(
        "https://raw.githubusercontent.com/anhph95/ifcbtools/main/analysis/community-variability/python/install_analysis_scripts.py"
    )
    .read()
    .decode("utf-8")
)
```

Open the copied scripts and build on them. Keep `ifcb_common.py` in the same
directory as the analysis scripts because each script imports it with
`from ifcb_common import ...`.

Each script also writes timestamped `.out.log` and `.err.log` files under
`<current-directory>/logs`. Use `--log-level DEBUG` for more detail or
`--log-dir PATH` to choose another log directory; terminal logging remains
enabled. Each log begins with the command-line arguments and resolved workflow
paths; secret-like parameter names are automatically redacted.
