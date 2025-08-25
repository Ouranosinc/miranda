from __future__ import annotations

from schema import And, Optional, Or, Regex, Schema

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


def _time_units_in_dimensions(time_dict: dict):
    """
    Check for time units in the dimensions.

    Parameters
    ----------
    time_dict : dict
        The Schema dictionary for the dimensions.

    Returns
    -------
    Schema
        The validated header Schema dictionary.

    Raises
    ------
    ValueError
        If the dimensions does not contain "time" as a key.
        If the time dimensions has neither "units" nor {"_units": true}
        If the time dimensions contains both "units" and {"_units": true}
    """
    if not time_dict or not isinstance(time_dict, dict):
        raise ValueError("'time' must be present and be a dictionary.")

    # Must contain either 'units' or '_units', but not both
    has_units = "units" in time_dict
    has_calculate_units = "_units" in time_dict

    if has_units and has_calculate_units:  # both true or both false
        raise ValueError("Time dimension may contain either 'units' or '_units', but not both")
    return time_dict


# Templates for different CF-compliant dimension types
# Time
cf_time_dimension_schema = And(
    Schema(
        {
            Optional("_cf_dimension_name"): "time",
            Optional("_ensure_correct_time"): Or(
                Schema(Regex(VALID_TIME_FREQUENCY_REGEX)),
                Schema({Regex(PROJECT_NAME_REGEX): Regex(VALID_TIME_FREQUENCY_REGEX)}),
            ),
            Optional("_strict_time"): bool,
            Optional("_units"): True,
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
            Optional("long_name"): Regex(r"^time$|^Time$"),
            "standard_name": "time",
            Optional("units"): Regex(TIME_UNITS_REGEX),
        },
        name="cf_time_dimension_schema",
    ),
    _time_units_in_dimensions,  # Validate time units
)


# Latitude
cf_lat_dimension_schema = Schema(
    {
        Optional("_cf_dimension_name"): "lat",
        Optional("_precision"): Or(
            int,
            Schema({Regex(PROJECT_NAME_REGEX): int}),
        ),
        "axis": "Y",
        Optional("bounds"): "lat_bnds",
        Optional("long_name"): Regex(r"^latitude$|^Latitude$"),
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
            Schema({Regex(PROJECT_NAME_REGEX): int}),
        ),
        "axis": "X",
        Optional("bounds"): "lon_bnds",
        Optional("long_name"): Regex(r"^longitude$|^Longitude$"),
        "standard_name": "longitude",
        "units": str,
    },
    name="cf_lon_dimension_schema",
)


# Dimensions
cf_dimensions_schema = Schema(
    {
        Optional("time"): cf_time_dimension_schema,
        Optional("lat"): cf_lat_dimension_schema,
        Optional("lon"): cf_lon_dimension_schema,
    },
    name="cf_dimensions_schema",
    ignore_extra_keys=True,
)
