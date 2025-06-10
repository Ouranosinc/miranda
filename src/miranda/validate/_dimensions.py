from __future__ import annotations

from schema import Optional

from .regex import (
    CELL_METHODS_REGEX,
    CF_CONVENTIONS_REGEX,
    PROJECT_NAME_REGEX,
    STANDARD_NAME_REGEX,
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
cf_time_dimension_schema = {
    Optional("_cf_dimension_name"): "time",
    Optional("_ensure_correct_time"): Or(
        Regex(VALID_TIME_FREQUENCIES),
        {Regex(PROJECT_NAME_REGEX): Regex(VALID_TIME_FREQUENCIES)},
    ),
    Optional("_strict_time"): bool,
    "axis": "T",
    "standard_name": "time",
    Optional("long_name"): Or("time", "Time"),
}

# Latitude
cf_lat_dimension_schema = {
    Optional("_cf_dimension_name"): "lat",
    Optional("_precision"): Or(
        int,
        {Regex(PROJECT_NAME_REGEX): int},
    ),
    "axis": "Y",
    "standard_name": "latitude",
    Optional("long_name"): Or("latitude", "Latitude"),
}

# Longitude
cf_lon_dimension_schema = {
    Optional("_cf_dimension_name"): "lon",
    Optional("_precision"): Or(
        int,
        {Regex(PROJECT_NAME_REGEX): int},
    ),
    "axis": "X",
    "standard_name": "longitude",
    Optional("long_name"): Or("longitude", "Longitude"),
}

# Dimensions
cf_dimensions_schema = {
    Optional("dimensions"): {
        Optional(str): Or(
            cf_time_dimension_schema,
            cf_lat_dimension_schema,
            cf_lon_dimension_schema,
        )
    }
}
