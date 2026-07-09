"""Taxonomy loading helpers."""

from __future__ import annotations

import argparse
import logging
import os
import re
from pathlib import Path
from typing import Callable, Sequence

import pandas as pd

LOGGER = logging.getLogger("ifcb")


def import_google_sheet(share_url: str, save_path: str | os.PathLike[str] | None = None) -> pd.DataFrame:
    """Import a Google Sheet as CSV and optionally save it locally."""
    id_match = re.search(r"/d/([a-zA-Z0-9-_]+)", share_url)
    if not id_match:
        raise ValueError("Spreadsheet ID not found in URL.")
    spreadsheet_id = id_match.group(1)

    gid_match = re.search(r"[?&]gid=(\d+)", share_url)
    gid = gid_match.group(1) if gid_match else "0"
    export_url = f"https://docs.google.com/spreadsheets/d/{spreadsheet_id}/export?format=csv&gid={gid}"

    LOGGER.info("Downloading taxonomy from Google Sheets: %s", export_url)
    df = pd.read_csv(export_url)
    taxonomic_levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
    present_levels = [col for col in taxonomic_levels if col in df.columns]
    if present_levels:
        df["Label"] = df[present_levels].ffill(axis=1).iloc[:, -1]

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(save_path, index=False)
        LOGGER.info("Saved taxonomy to: %s", save_path)
    return df


DEFAULT_TAXONOMIC_LEVELS = ("Phylum", "Class", "Order", "Family", "Genus", "Species")


def add_taxonomy_label(
    taxonomy: pd.DataFrame,
    taxonomic_levels: Sequence[str] = DEFAULT_TAXONOMIC_LEVELS,
    label_col: str = "Label",
) -> pd.DataFrame:
    """Create a filled taxon label from the deepest available taxonomic level.

    For each annotation j, the label is the last non-missing value along:
    Phylum -> Class -> Order -> Family -> Genus -> Species.
    """
    taxonomy = taxonomy.copy()
    for col in taxonomy.columns:
        if pd.api.types.is_object_dtype(taxonomy[col]) or pd.api.types.is_string_dtype(taxonomy[col]):
            taxonomy[col] = taxonomy[col].map(lambda x: x.strip() if isinstance(x, str) else x)
    levels = [col for col in taxonomic_levels if col in taxonomy.columns]
    if not levels:
        raise ValueError("taxonomy must include at least one taxonomic level column.")
    taxonomy[label_col] = taxonomy[levels].ffill(axis=1).iloc[:, -1]
    return taxonomy


def map_taxa_to_label(
    df: pd.DataFrame,
    taxonomy: pd.DataFrame,
    from_level: str = "Annotations",
    to_level: str = "Label",
    aggfunc: str = "sum",
) -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    """Aggregate annotation columns to a taxonomic label while keeping metadata.

    Mathematically, if several annotation columns j map to the same label g,
    the filled table stores X_g = sum_j X_j for each sample row.
    """
    taxonomy = taxonomy.copy()
    for col in taxonomy.columns:
        if pd.api.types.is_object_dtype(taxonomy[col]) or pd.api.types.is_string_dtype(taxonomy[col]):
            taxonomy[col] = taxonomy[col].map(lambda x: x.strip() if isinstance(x, str) else x)
    if to_level not in taxonomy.columns:
        levels = [col for col in DEFAULT_TAXONOMIC_LEVELS if col in taxonomy.columns]
        if not levels:
            raise ValueError("taxonomy must include at least one taxonomic level column.")
        taxonomy[to_level] = taxonomy[levels].ffill(axis=1).iloc[:, -1]
    if from_level not in taxonomy.columns or to_level not in taxonomy.columns:
        raise ValueError(f"taxonomy must include '{from_level}' and '{to_level}'.")

    mapping = taxonomy.set_index(from_level)[to_level].dropna().to_dict()
    mapped_taxa = [col for col in df.columns if col in mapping]
    metadata_cols = [col for col in df.columns if col not in mapped_taxa]

    taxa = df.loc[:, mapped_taxa].apply(pd.to_numeric, errors="coerce").copy()
    taxa.columns = [mapping[col] for col in mapped_taxa]
    taxa_agg = taxa.T.groupby(level=0).agg(aggfunc).T

    out = pd.concat([df.loc[:, metadata_cols].copy(), taxa_agg], axis=1)
    LOGGER.info("Mapped %s annotation columns to %s taxon labels", len(mapped_taxa), len(taxa_agg.columns))
    return out, sorted(taxa_agg.columns.tolist()), taxonomy


def taxon_mapping(
    df_full: pd.DataFrame,
    taxonomy: pd.DataFrame | str | Path = "taxonomy.csv",
    from_level: str = "Annotations",
    to_level: str = "Genus",
    aggfunc: str | Callable = "sum",
) -> tuple[pd.DataFrame, list[str]]:
    """Aggregate taxon columns from one taxonomy level to another.

    Parameters
    ----------
    df_full
        Input data with metadata columns and taxon columns named at
        ``from_level``.
    taxonomy
        Taxonomy table, or a path to a CSV taxonomy table. The table must
        include ``from_level`` and ``to_level`` columns.
    from_level
        Taxonomy level currently used by taxon columns in ``df_full``.
    to_level
        Taxonomy level to aggregate taxon columns to.
    aggfunc
        Aggregation passed to :meth:`pandas.core.groupby.DataFrameGroupBy.agg`,
        such as ``"sum"`` or ``"mean"``.

    Returns
    -------
    tuple[pandas.DataFrame, list[str]]
        A dataframe containing metadata columns plus aggregated taxon columns,
        and the sorted list of aggregated taxon column names.
    """
    if not isinstance(taxonomy, pd.DataFrame):
        taxonomy = pd.read_csv(taxonomy)
    else:
        taxonomy = taxonomy.copy()

    if from_level not in taxonomy.columns or to_level not in taxonomy.columns:
        raise ValueError(
            f"Taxonomy table must include columns {from_level!r} and {to_level!r}"
        )

    taxonomy = taxonomy.apply(
        lambda col: col.map(lambda value: value.strip() if isinstance(value, str) else value)
    )

    mapping = taxonomy.set_index(from_level)[to_level].dropna().to_dict()
    df_columns_stripped = {
        col: col.strip() if isinstance(col, str) else col for col in df_full.columns
    }
    mapped_taxa = [
        col for col, stripped_col in df_columns_stripped.items() if stripped_col in mapping
    ]
    unmapped_cols = [col for col in df_full.columns if col not in mapped_taxa]

    df_mapped = df_full[mapped_taxa].copy()
    df_mapped.columns = [mapping[df_columns_stripped[col]] for col in mapped_taxa]
    df_agg = df_mapped.T.groupby(level=0).agg(aggfunc).T

    df_result = pd.concat([df_full[unmapped_cols], df_agg], axis=1)
    return df_result, sorted(df_agg.columns.tolist())


def taxonomy_mapping(
    df_full: pd.DataFrame,
    taxonomy: pd.DataFrame | str | Path = "taxonomy.csv",
    from_level: str = "Annotations",
    to_level: str = "Genus",
    aggfunc: str | Callable = "sum",
) -> tuple[pd.DataFrame, list[str]]:
    """Alias for :func:`taxon_mapping` using the more explicit function name."""
    return taxon_mapping(
        df_full,
        taxonomy=taxonomy,
        from_level=from_level,
        to_level=to_level,
        aggfunc=aggfunc,
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for taxonomy mapping."""
    parser = argparse.ArgumentParser(
        description="Aggregate IFCB taxon columns from one taxonomy level to another.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", dest="input_file", required=True, help="Input CSV file.")
    parser.add_argument("--taxonomy-file", default="taxonomy.csv", help="Taxonomy CSV file.")
    parser.add_argument("--output", "--output-file", dest="output_file", required=True, help="Output CSV file.")
    parser.add_argument("--from-level", default="Annotations", help="Taxonomy level used by input taxon columns.")
    parser.add_argument("--to-level", default="Genus", help="Taxonomy level to aggregate taxon columns to.")
    parser.add_argument("--aggfunc", default="sum", help="Pandas groupby aggregation function.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run taxonomy mapping from the command line."""
    args = parse_args(argv)
    input_path = Path(args.input_file)
    taxonomy_path = Path(args.taxonomy_file)
    output_path = Path(args.output_file)

    try:
        if not input_path.exists():
            raise FileNotFoundError(f"Input file does not exist: {input_path}")
        if not taxonomy_path.exists():
            raise FileNotFoundError(f"Taxonomy file does not exist: {taxonomy_path}")

        df = pd.read_csv(input_path, low_memory=False)
        mapped, taxa = taxon_mapping(
            df,
            taxonomy=taxonomy_path,
            from_level=args.from_level,
            to_level=args.to_level,
            aggfunc=args.aggfunc,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        mapped.to_csv(output_path, index=False)
    except Exception:
        LOGGER.exception("IFCB taxonomy mapping failed")
        return 1

    LOGGER.info("Mapped %s taxon columns. Output: %s", len(taxa), output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
