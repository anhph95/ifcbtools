"""Main IFCB processing orchestrator."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys
from typing import Sequence

import pandas as pd

from ifcb.logging import log_run_configuration, redact_command_line, setup_logging

from . import add
from .clean import filter_and_normalize
from .taxonomy import import_google_sheet

LOGGER = logging.getLogger("ifcb")


def default_output_path(
    input_file: str | Path,
    *,
    clean: bool = False,
    station: bool = False,
    bottle: bool = False,
    nutrient: bool = False,
) -> Path:
    """Append selected operation suffixes to an input filename."""
    input_path = Path(input_file)
    suffixes = []
    if clean:
        suffixes.append("clean")
    if station:
        suffixes.append("station")
    if bottle:
        suffixes.append("bottle")
    if nutrient:
        suffixes.append("nutrient")
    return input_path.with_name(f"{input_path.stem}{''.join(f'_{suffix}' for suffix in suffixes)}{input_path.suffix}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process one IFCB CSV through selected cleaning and enrichment steps.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", dest="input_file", required=True, help="Input CSV file.")
    parser.add_argument("--output", "--output-file", dest="output_file", default=None, help="Final output CSV file.")
    parser.add_argument("--all", action="store_true", help="Run clean, station, bottle, and nutrient unless explicitly disabled.")
    parser.add_argument("--clean", action=argparse.BooleanOptionalAction, default=None, help="Quality-control metadata, aggregate casts, and normalize values.")
    parser.add_argument("--station", action=argparse.BooleanOptionalAction, default=None, help="Assign nearest NES-LTER station fields.")
    parser.add_argument("--bottle", action=argparse.BooleanOptionalAction, default=None, help="Merge CTD bottle fields.")
    parser.add_argument("--nutrient", action=argparse.BooleanOptionalAction, default=None, help="Merge nutrient fields.")
    parser.add_argument(
        "--taxonomy-file",
        default=argparse.SUPPRESS,
        help="Taxonomy CSV. Omit to use ifcb_taxonomy.csv beside the input; missing taxonomy is downloaded during --clean.",
    )
    parser.add_argument(
        "--data-type",
        choices=["count", "carbon"],
        default=argparse.SUPPRESS,
        help="Product type for --clean normalization. Omit to infer from input filename.",
    )
    parser.add_argument(
        "--taxonomy-url",
        default=argparse.SUPPRESS,
        help="Google Sheet taxonomy URL used only when the taxonomy CSV is missing.",
    )
    parser.add_argument(
        "--station-reference",
        default=argparse.SUPPRESS,
        help="Station reference CSV or URL. Omit to use the station-assignment function default.",
    )
    parser.add_argument(
        "--max-station-distance-km",
        type=float,
        default=argparse.SUPPRESS,
        help="Maximum station assignment distance in km. Omit to use the station-assignment function default; use a large value for effectively no limit.",
    )
    parser.add_argument(
        "--bottle-url-template",
        default=argparse.SUPPRESS,
        help="CTD bottle CSV URL template. Omit to use the bottle-merge function default.",
    )
    parser.add_argument(
        "--nutrient-source",
        default=argparse.SUPPRESS,
        help="Nutrient CSV path or URL. Omit to use the nutrient-merge function default.",
    )
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO", help="Logging level. Default: INFO.")
    parser.add_argument("--log-dir", default=None, help="Directory for timestamped logs. Omit to use ./logs.")
    args = parser.parse_args(argv)

    if args.all:
        for name in ("clean", "station", "bottle", "nutrient"):
            if getattr(args, name) is None:
                setattr(args, name, True)
    else:
        for name in ("clean", "station", "bottle", "nutrient"):
            if getattr(args, name) is None:
                setattr(args, name, False)

    if not any((args.clean, args.station, args.bottle, args.nutrient)):
        parser.error("select at least one operation with --all or explicit True flags.")
    return args


def process(
    input_file: str | Path,
    output_file: str | Path | None = None,
    *,
    clean: bool = False,
    station: bool = False,
    bottle: bool = False,
    nutrient: bool = False,
    taxonomy_file: str | Path | None = None,
    data_type: str | None = None,
    taxonomy_url: str = (
        "https://docs.google.com/spreadsheets/d/"
        "1mTmL3VN3TlNvkz4uJpDp34QQb5ImUlC4J7FpxKs4nFI/edit?gid=1521292620"
    ),
    station_reference: str | Path | None = "https://nes-lter-api.whoi.edu/api/stations/file.csv",
    max_station_distance_km: float | None = 2.0,
    bottle_url_template: str = "https://nes-lter-api.whoi.edu/api/ctd/bottles/{cruise}.csv",
    nutrient_source: str | Path = "https://nes-lter-api.whoi.edu/api/nut/all.csv",
) -> Path:
    """Run selected IFCB processing steps in pipeline order."""
    input_path = Path(input_file)
    final_path = Path(output_file) if output_file is not None else default_output_path(
        input_path,
        clean=clean,
        station=station,
        bottle=bottle,
        nutrient=nutrient,
    )
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")
    df = pd.read_csv(input_path, low_memory=False)
    taxonomy = None

    if clean:
        taxonomy_path = Path(taxonomy_file) if taxonomy_file is not None else input_path.parent / "ifcb_taxonomy.csv"
        if not taxonomy_path.exists():
            taxonomy = import_google_sheet(taxonomy_url, save_path=taxonomy_path)

        missing = [path for path in (taxonomy_path,) if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Selected input file(s) do not exist: {missing}")

        if taxonomy is None:
            taxonomy = pd.read_csv(taxonomy_path, low_memory=False)

        selected_data_type = data_type
        if selected_data_type is None:
            stem = input_path.stem.lower()
            if "count" in stem:
                selected_data_type = "count"
            elif "carbon" in stem:
                selected_data_type = "carbon"
            else:
                raise ValueError(
                    "Cannot infer IFCB data type from input filename. "
                    "Re-run with --data-type count or --data-type carbon."
                )
        if selected_data_type == "count":
            scaling_factor = 1000.0
        elif selected_data_type == "carbon":
            scaling_factor = 0.001
        else:
            raise ValueError("data_type must be one of: count, carbon")

        df = filter_and_normalize(df, taxonomy, scaling_factor=scaling_factor)

    if station:
        df = add.nearest_station(
            df,
            station_reference=station_reference,
            max_station_distance_km=max_station_distance_km,
        )
    if bottle:
        df = add.bottle(df, bottle_url_template=bottle_url_template)
    if nutrient:
        df = add.nutrient(df, nutrient_source=nutrient_source)

    final_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(final_path, index=False)
    LOGGER.info("Saved processed data to: %s", final_path)
    return final_path


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    log_dir = Path(args.log_dir) if args.log_dir is not None else Path.cwd() / "logs"
    setup_logging(log_dir=log_dir, name="ifcb_process", level=getattr(logging, args.log_level))
    process_options = {
        name: getattr(args, name)
        for name in (
            "taxonomy_file",
            "data_type",
            "taxonomy_url",
            "station_reference",
            "max_station_distance_km",
            "bottle_url_template",
            "nutrient_source",
        )
        if hasattr(args, name)
    }
    log_run_configuration(
        LOGGER,
        {
            "command": redact_command_line(sys.argv),
            "input_file": Path(args.input_file).resolve(),
            "output_file": Path(args.output_file).resolve() if args.output_file is not None else None,
            "clean": args.clean,
            "station": args.station,
            "bottle": args.bottle,
            "nutrient": args.nutrient,
            "all": args.all,
            **process_options,
        },
    )
    try:
        output = process(
            input_file=args.input_file,
            output_file=args.output_file,
            clean=args.clean,
            station=args.station,
            bottle=args.bottle,
            nutrient=args.nutrient,
            **process_options,
        )
    except Exception:
        LOGGER.exception("IFCB process failed")
        return 1
    LOGGER.info("IFCB process completed. Output: %s", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
