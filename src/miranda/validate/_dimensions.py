from __future__ import annotations

from schema import Optional, Or, Regex, Schema

from ._regex import (
    PROJECT_NAME_REGEX,
    TIME_UNITS_REGEX,
    VALID_TIME_FREQUENCY_REGEX,
)

__all__ = ["cf_dimensions_schema"]

CF_DIMENSION_TREATMENTS = [
    "_cf_dimension_name",
    "_ensure_correct_time_frequency",
    "_offset_time",
    "_precision",
    "_strict_time",
]


# Templates for different CF-compliant dimension types
# Time
cf_time_dimension_schema = Schema(
    {
        Optional("_cf_dimension_name"): "time",
        Optional("_ensure_correct_time"): Or(
            Regex(VALID_TIME_FREQUENCY_REGEX),
            {Regex(PROJECT_NAME_REGEX): Regex(VALID_TIME_FREQUENCY_REGEX)},
        ),
        Optional("_strict_time"): bool,
        "axis": "T",
        Optional("bounds"): "time_bnds",
        Optional("calendar"): Or(
            "standard",
            "gregorian",
            "proleptic_gregorian",
            "noleap",
            "365_day",
            "360_day",
        ),
        Optional("long_name"): Or("time", "Time"),
        "standard_name": "time",
        "units": Regex(TIME_UNITS_REGEX),
    },
    name="cf_time_dimension_schema",
)


# Latitude
cf_lat_dimension_schema = Schema(
    {
        Optional("_cf_dimension_name"): "lat",
        Optional("_precision"): Or(
            int,
            {Regex(PROJECT_NAME_REGEX): int},
        ),
        "axis": "Y",
        Optional("bounds"): "lat_bnds",
        Optional("long_name"): Or("latitude", "Latitude"),
        "standard_name": "latitude",
        "units": str,
    },
    name="cf_lat_dimension_schema",
)


# Longitude
cf_lon_dimension_schema = Schema(
    {
        Optional("_cf_dimension_name"): "lon",
        Optional("_precision"): Or(
            int,
            {Regex(PROJECT_NAME_REGEX): int},
        ),
        "axis": "X",
        Optional("bounds"): "lon_bnds",
        Optional("long_name"): Or("longitude", "Longitude"),
        "standard_name": "longitude",
        "units": str,
    },
    name="cf_lon_dimension_schema",
)


# Dimensions
cf_dimensions_schema = Schema(
    {
        Optional("dimensions"): {
            Optional("time"): cf_time_dimension_schema,
            Optional("lat"): cf_lat_dimension_schema,
            Optional("lon"): cf_lon_dimension_schema,
        }
    },
    name="cf_dimensions_schema",
)
