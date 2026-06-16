"""Command-line interface for IFCB data processing."""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Sequence

from .constants import DEFAULT_BOTTLE_URL_TEMPLATE, DEFAULT_DATASET, DEFAULT_TAXONOMY_URL
from .fill import DEFAULT_NUTRIENT_URL
from .process import process

LOGGER = logging.getLogger("ifcb.process.neslter")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process NES-LTER IFCB CSV files exported by MATLAB.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_data_path",
        nargs="?",
        default=None,
        help="Directory containing IFCB input CSV files. Defaults to data/<dataset>.",
    )
    parser.add_argument("--dataset", default=DEFAULT_DATASET, help="Dataset folder under repo-local data/.")
    parser.add_argument("-o", "--output-dir", default=None, help="Directory for cleaned output CSVs.")
    parser.add_argument("--sample-type", nargs="+", default=None, help="sample_type values to keep.")
    parser.add_argument(
        "--data-type",
        choices=["count", "carbon"],
        nargs="+",
        default=["count", "carbon"],
        help="Raw data types to process.",
    )
    parser.add_argument(
        "--download-taxonomy-if-missing",
        action="store_true",
        default=True,
        help="Download ifcb_taxonomy.csv from Google Sheets when missing.",
    )
    parser.add_argument(
        "--no-download-taxonomy",
        action="store_false",
        dest="download_taxonomy_if_missing",
        help="Require an existing ifcb_taxonomy.csv instead of downloading it.",
    )
    parser.add_argument("--taxonomy-url", default=DEFAULT_TAXONOMY_URL)
    parser.add_argument("--station-reference", default=None, help="Station reference CSV path.")
    parser.add_argument("--max-station-distance-km", type=float, default=2.0)
    parser.add_argument("--no-station-distance-limit", action="store_true")
    parser.add_argument("--bottle-url-template", default=DEFAULT_BOTTLE_URL_TEMPLATE)
    parser.add_argument("--nutrient-source", default=DEFAULT_NUTRIENT_URL)
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    max_distance = None if args.no_station_distance_limit else args.max_station_distance_km

    try:
        outputs = process(
            input_dir=args.input_data_path,
            output_dir=args.output_dir,
            dataset=args.dataset,
            sample_type=args.sample_type,
            download_taxonomy_if_missing=args.download_taxonomy_if_missing,
            taxonomy_url=args.taxonomy_url,
            station_reference=args.station_reference,
            max_station_distance_km=max_distance,
            bottle_url_template=args.bottle_url_template,
            nutrient_source=args.nutrient_source,
            data_types=args.data_type,
        )
    except Exception as exc:
        LOGGER.error("Process failed: %s", exc)
        return 1

    LOGGER.info("Process completed. Outputs: %s", ", ".join(str(path) for path in outputs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
