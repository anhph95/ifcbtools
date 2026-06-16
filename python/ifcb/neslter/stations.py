"""Station lookup utilities.

This module is intentionally independent from the rest of the IFCB process so it
can be imported and tested on its own.
"""

from __future__ import annotations

import os
from typing import Sequence

import numpy as np
import pandas as pd
from geopy.distance import geodesic as geo_distance

from .constants import DEFAULT_STATION_REF_URL


def load_station_reference(station_reference: pd.DataFrame | str | os.PathLike[str] | None = None) -> pd.DataFrame:
    """Load and normalize the NES-LTER station reference table."""
    if station_reference is None:
        st = pd.read_csv(DEFAULT_STATION_REF_URL)
    elif isinstance(station_reference, (str, os.PathLike)):
        st = pd.read_csv(station_reference)
    else:
        st = station_reference.copy()

    st.columns = st.columns.str.strip()
    required = ["station", "startDate", "endDate", "decimalLatitude", "decimalLongitude"]
    missing = [col for col in required if col not in st.columns]
    if missing:
        raise ValueError(f"station_reference missing required columns: {missing}")

    st["station"] = st["station"].astype(str).str.strip()
    st["startDate"] = pd.to_datetime(st["startDate"], errors="coerce", utc=True).dt.tz_localize(None)
    st["endDate"] = st["endDate"].replace("current", pd.NA)
    st["endDate"] = pd.to_datetime(st["endDate"], errors="coerce", utc=True).dt.tz_localize(None)
    st["endDate_filled"] = st["endDate"].fillna(pd.Timestamp("2100-12-31"))
    st["decimalLatitude"] = pd.to_numeric(st["decimalLatitude"], errors="coerce")
    st["decimalLongitude"] = pd.to_numeric(st["decimalLongitude"], errors="coerce")
    return st.dropna(subset=["station", "startDate", "decimalLatitude", "decimalLongitude"]).copy()


def active_stations(station_reference: pd.DataFrame, timestamp: object) -> pd.DataFrame:
    """Return stations active at a timestamp."""
    ts = pd.to_datetime(timestamp, errors="coerce", utc=True)
    if pd.isna(ts):
        raise ValueError("timestamp is missing or invalid")
    ts = ts.tz_localize(None)

    active = station_reference[
        (station_reference["startDate"] <= ts) & (ts <= station_reference["endDate_filled"])
    ].copy()
    if active.empty:
        raise ValueError(f"No active stations found for timestamp {ts}")
    return active


def station_distances(
    lat: float,
    lon: float,
    timestamp: object,
    station_reference: pd.DataFrame | str | os.PathLike[str] | None = None,
) -> pd.Series:
    """Calculate distance in km from one coordinate to active stations."""
    st = load_station_reference(station_reference)
    active = active_stations(st, timestamp)
    distances = []
    index = []
    for row in active.itertuples():
        distances.append(geo_distance((lat, lon), (row.decimalLatitude, row.decimalLongitude)).km)
        index.append(row.Index)
    return pd.Series(distances, index=index)


def nearest_station(
    lat: float,
    lon: float,
    timestamp: object,
    station_reference: pd.DataFrame | str | os.PathLike[str] | None = None,
    max_distance_km: float | None = 2.0,
) -> tuple[object, float | object]:
    """Return the nearest active station name and distance in km.

    This is the standalone helper function to use outside the full data process.
    """
    st = load_station_reference(station_reference)
    distances = station_distances(lat=lat, lon=lon, timestamp=timestamp, station_reference=st)
    idx = distances.idxmin()
    d_km = float(distances.loc[idx])
    if max_distance_km is not None and d_km > max_distance_km:
        return np.nan, np.nan
    return st.loc[idx, "station"], d_km


def assign_nearest_stations(
    df: pd.DataFrame,
    station_reference: pd.DataFrame | str | os.PathLike[str] | None = None,
    lat_col: str = "latitude",
    lon_col: str = "longitude",
    time_col: str = "sample_time",
    max_distance_km: float | None = 2.0,
    output_station_col: str = "nearest_station",
    output_distance_col: str = "station_distance",
) -> pd.DataFrame:
    """Assign nearest station columns to a DataFrame."""
    st = load_station_reference(station_reference)
    out = df.copy()
    names: list[object] = []
    distances: list[object] = []

    for row in out.itertuples(index=False):
        lat = getattr(row, lat_col)
        lon = getattr(row, lon_col)
        ts = getattr(row, time_col)
        if pd.isna(lat) or pd.isna(lon) or pd.isna(ts):
            names.append(np.nan)
            distances.append(np.nan)
            continue
        name, dist = nearest_station(lat, lon, ts, station_reference=st, max_distance_km=max_distance_km)
        names.append(name)
        distances.append(dist)

    out[output_station_col] = names
    out[output_distance_col] = distances
    return out


class StationLocator:
    """Small convenience wrapper around the standalone station functions."""

    def __init__(
        self,
        station_reference: pd.DataFrame | str | os.PathLike[str] | None = None,
        max_distance_km: float | None = 2.0,
    ) -> None:
        self.station_reference = load_station_reference(station_reference)
        self.max_distance_km = max_distance_km

    def active_stations(self, timestamp: object) -> pd.DataFrame:
        return active_stations(self.station_reference, timestamp)

    def station_distances(self, lat: float, lon: float, timestamp: object) -> pd.Series:
        return station_distances(lat, lon, timestamp, self.station_reference)

    def nearest_station(
        self,
        lat: float,
        lon: float,
        timestamp: object,
        max_distance_km: float | None = None,
    ) -> tuple[object, float | object]:
        threshold = self.max_distance_km if max_distance_km is None else max_distance_km
        return nearest_station(lat, lon, timestamp, self.station_reference, threshold)

    def nearest_stations(
        self,
        df: pd.DataFrame,
        lat_col: str = "latitude",
        lon_col: str = "longitude",
        time_col: str = "sample_time",
        max_distance_km: float | None = None,
    ) -> tuple[list[object], list[object]]:
        threshold = self.max_distance_km if max_distance_km is None else max_distance_km
        assigned = assign_nearest_stations(
            df,
            station_reference=self.station_reference,
            lat_col=lat_col,
            lon_col=lon_col,
            time_col=time_col,
            max_distance_km=threshold,
        )
        return assigned["nearest_station"].tolist(), assigned["station_distance"].tolist()
