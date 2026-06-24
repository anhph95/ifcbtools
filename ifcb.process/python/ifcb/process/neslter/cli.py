"""Command-line interface for IFCB data processing."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys
from typing import Sequence

from ifcb.process.logging_utils import log_run_configuration, redact_command_line, setup_logging

from .constants import DEFAULT_BOTTLE_URL_TEMPLATE, DEFAULT_TAXONOMY_URL
from .fill import DEFAULT_NUTRIENT_URL
from .process import default_output_path, process

LOGGER = logging.getLogger("ifcb.process.neslter")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process NES-LTER IFCB CSV files exported by MATLAB.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_file",
        help="Input CSV file to process.",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        default=None,
        help="Output CSV file. By default, operation suffixes are added to the input filename.",
    )
    parser.add_argument(
        "--metadata-file",
        default=None,
        help="Metadata CSV used by --clean. Defaults to ifcb_metadata.csv beside the input file.",
    )
    parser.add_argument(
        "--taxonomy-file",
        default=None,
        help="Taxonomy CSV used by --clean. Defaults to ifcb_taxonomy.csv beside the input file.",
    )
    operations = parser.add_argument_group("processing operations")
    operations.add_argument(
        "--all",
        action="store_true",
        help="Run clean, station assignment, bottle merge, and nutrient merge.",
    )
    operations.add_argument(
        "--clean",
        action="store_true",
        help="Clean raw IFCB data and create the output CSV.",
    )
    operations.add_argument(
        "--add-station",
        action="store_true",
        help="Add nearest-station fields to the output CSV.",
    )
    operations.add_argument(
        "--merge-bottle",
        action="store_true",
        help="Merge CTD bottle fields into the output CSV.",
    )
    operations.add_argument(
        "--merge-nutrient",
        action="store_true",
        help="Merge nutrient fields into the output CSV.",
    )
    parser.add_argument("--sample-type", nargs="+", default=None, help="sample_type values to keep.")
    parser.add_argument(
        "--download-taxonomy-if-missing",
        action="store_true",
        default=True,
        help="Download the selected taxonomy file from Google Sheets when missing.",
    )
    parser.add_argument(
        "--no-download-taxonomy",
        action="store_false",
        dest="download_taxonomy_if_missing",
        help="Require the selected taxonomy file instead of downloading it.",
    )
    parser.add_argument("--taxonomy-url", default=DEFAULT_TAXONOMY_URL)
    parser.add_argument("--station-reference", default=None, help="Station reference CSV path.")
    parser.add_argument("--max-station-distance-km", type=float, default=2.0)
    parser.add_argument("--no-station-distance-limit", action="store_true")
    parser.add_argument("--bottle-url-template", default=DEFAULT_BOTTLE_URL_TEMPLATE)
    parser.add_argument("--nutrient-source", default=DEFAULT_NUTRIENT_URL)
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    parser.add_argument("--log-dir", default=None, help="Directory for timestamped workflow log files.")
    args = parser.parse_args(argv)
    if not any((args.all, args.clean, args.add_station, args.merge_bottle, args.merge_nutrient)):
        parser.error("select at least one operation: --all, --clean, --add-station, --merge-bottle, or --merge-nutrient")
    if args.all:
        args.clean = True
        args.add_station = True
        args.merge_bottle = True
        args.merge_nutrient = True
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input_file)
    output_path = (
        Path(args.output_file)
        if args.output_file is not None
        else default_output_path(
            input_path,
            clean=args.clean,
            add_station=args.add_station,
            merge_bottle_data=args.merge_bottle,
            merge_nutrient=args.merge_nutrient,
        )
    )
    metadata_path = Path(args.metadata_file) if args.metadata_file is not None else input_path.parent / "ifcb_metadata.csv"
    taxonomy_path = Path(args.taxonomy_file) if args.taxonomy_file is not None else input_path.parent / "ifcb_taxonomy.csv"
    log_dir = Path(args.log_dir) if args.log_dir is not None else Path.cwd() / "logs"
    setup_logging(log_dir=log_dir, name="ifcb_process", level=getattr(logging, args.log_level))
    max_distance = None if args.no_station_distance_limit else args.max_station_distance_km
    LOGGER.info("Starting IFCB processing for %s", input_path)
    log_run_configuration(
        LOGGER,
        {
            "command": redact_command_line(sys.argv),
            "input_file": input_path.resolve(),
            "output_file": output_path.resolve(),
            "sample_type": args.sample_type,
            "all": args.all,
            "metadata_file": metadata_path.resolve(),
            "taxonomy_file": taxonomy_path.resolve(),
            "clean": args.clean,
            "add_station": args.add_station,
            "merge_bottle": args.merge_bottle,
            "merge_nutrient": args.merge_nutrient,
            "download_taxonomy_if_missing": args.download_taxonomy_if_missing,
            "taxonomy_url": args.taxonomy_url,
            "station_reference": args.station_reference,
            "max_station_distance_km": max_distance,
            "bottle_url_template": args.bottle_url_template,
            "nutrient_source": args.nutrient_source,
            "log_level": args.log_level,
            "log_dir": log_dir.resolve(),
        },
    )

    try:
        output = process(
            input_file=input_path,
            output_file=output_path,
            sample_type=args.sample_type,
            download_taxonomy_if_missing=args.download_taxonomy_if_missing,
            taxonomy_url=args.taxonomy_url,
            station_reference=args.station_reference,
            max_station_distance_km=max_distance,
            bottle_url_template=args.bottle_url_template,
            nutrient_source=args.nutrient_source,
            metadata_file=metadata_path,
            taxonomy_file=taxonomy_path,
            clean=args.clean,
            add_station=args.add_station,
            merge_bottle_data=args.merge_bottle,
            merge_nutrient=args.merge_nutrient,
        )
    except Exception:
        LOGGER.exception("Process failed")
        return 1

    LOGGER.info("Process completed. Output: %s", output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
