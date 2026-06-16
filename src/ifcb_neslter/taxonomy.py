"""Taxonomy loading helpers."""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path

import pandas as pd

LOGGER = logging.getLogger("ifcb_neslter")


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
        df["Mix_Taxa"] = df[present_levels].ffill(axis=1).iloc[:, -1]

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(save_path, index=False)
        LOGGER.info("Saved taxonomy to: %s", save_path)
    return df
