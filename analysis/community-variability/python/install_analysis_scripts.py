"""Download editable IFCB Python analysis starter scripts into the cwd."""

from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen


def install_ifcb_python_analysis_scripts(destdir: str | Path = ".", overwrite: bool = False) -> list[str]:
    base_url = (
        "https://raw.githubusercontent.com/anhph95/ifcbtools/main/"
        "analysis/community-variability/python/scripts"
    )
    script_files = [
        "ifcb_common.py",
        "ifcb_single_season.py",
        "ifcb_power_analysis.py",
        "ifcb_sensitivity_analysis.py",
        "ifcb_seasonal_comparison.py",
    ]

    target = Path(destdir).expanduser().resolve()
    print("Installing Python analysis scripts into:")
    print(f"  {target}")

    for script_file in script_files:
        destfile = target / script_file
        if destfile.exists() and not overwrite:
            print(f"Skipping existing file: {script_file}")
            continue

        url = f"{base_url}/{script_file}"
        with urlopen(url) as response:
            contents = response.read()
        destfile.write_bytes(contents)
        print(f"Downloaded: {script_file}")

    print("Done. Open the copied scripts in this working directory and build on them.")
    return script_files


install_ifcb_python_analysis_scripts(overwrite=bool(globals().get("IFCB_ANALYSIS_OVERWRITE", False)))
