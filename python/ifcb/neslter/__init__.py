"""Tools for processing NES-LTER IFCB data."""

from .community_variability import (
    BD_alpha,
    BD_gamma,
    BD_phi,
    BD_spatial_weighted,
    COMMUNITY_METRIC_ORDER,
    CV_alpha,
    CV_gamma,
    CV_phi,
    CommunityArray,
    add_baseline_delta,
    bootstrap_by_dimension,
    calc_metacommunity_metrics,
    calc_spatial_bd_by_time,
    leave_one_out,
    make_community_array,
    spatial_bd_by_time,
    wide_metric_table,
)
from .taxonomy import import_google_sheet

from .fill import make_filled_dataset

try:
    from .process import matlab_export_data_dir, process, process_data_type
    from .stations import StationLocator, assign_nearest_stations, nearest_station
except ModuleNotFoundError:
    matlab_export_data_dir = None
    process = None
    process_data_type = None
    StationLocator = None
    assign_nearest_stations = None
    nearest_station = None

__all__ = [
    "BD_alpha",
    "BD_gamma",
    "BD_phi",
    "BD_spatial_weighted",
    "COMMUNITY_METRIC_ORDER",
    "CV_alpha",
    "CV_gamma",
    "CV_phi",
    "CommunityArray",
    "StationLocator",
    "add_baseline_delta",
    "assign_nearest_stations",
    "bootstrap_by_dimension",
    "calc_metacommunity_metrics",
    "calc_spatial_bd_by_time",
    "import_google_sheet",
    "leave_one_out",
    "matlab_export_data_dir",
    "make_filled_dataset",
    "make_community_array",
    "nearest_station",
    "process",
    "process_data_type",
    "spatial_bd_by_time",
    "wide_metric_table",
]
