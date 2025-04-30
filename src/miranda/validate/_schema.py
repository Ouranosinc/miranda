from __future__ import annotations

import json
import logging
from pathlib import Path

from pydantic_core.core_schema import DictSchema
from schema import And, Optional, Or, Regex, Schema, SchemaError

from ._dimensions import cf_dimensions_schema
from .utils import (
    CELL_METHODS_REGEX,
    CF_CONVENTIONS_REGEX,
    PROJECT_NAME_REGEX,
    STANDARD_NAME_REGEX,
    VALID_TIME_FREQUENCIES,
)

__all__ = ["validate_json"]


DIMENSION_TREATMENTS = ["_cf_dimension_name", "", ""]

VARIABLE_TREATMENTS = [
    "_corrected_units",
    "_cf_units_conversion",
    "_clip_values",
    "_correct_unit_names",
    "_invert_value_sign",
    "_transform_values",
    "_units_context",
    "_variable_conversion",
]


def _institution_in_header(header_dict):
    """
    Check for institution metadata in header.

    Parameters
    ----------
    header_dict : Schema
        The Schema dictionary for the header.

    Returns
    -------
    Schema
        The validated header Schema dictionary.

    Raises
    ------
    ValueError
        If the Header has neither "institution" or {"_map_attrs": {str : "institution"}
        If the Header contains both "institution" and {"_map_attrs": {str : "institution"}
    """
    # Must contain either 'institution' or '_map_attrs', but not both
    has_institution = "institution" in header_dict
    has_map_attrs = "_map_attrs" in header_dict

    if has_institution == has_map_attrs:  # both true or both false
        raise ValueError(
            "Header must contain either 'institution' or '_map_attrs', but not both"
        )

    if has_institution:
        if not isinstance(header_dict["institution"], str):
            raise ValueError("'institution' must be a string")

    if has_map_attrs:
        map_attrs = header_dict["_map_attrs"]
        if not isinstance(map_attrs, dict) or not any(
            isinstance(k, str) and v == "institution" for k, v in map_attrs.items()
        ):
            raise ValueError("'_map_attrs' must be a dict of {str: 'institution'}")

    return header_dict


# Template for
header_schema = {
    "Header": And(
        {
            "Conventions": Regex(CF_CONVENTIONS_REGEX),
            "source": str,
            "type": str,
            "processing_level": Or("raw", "biasadjusted"),
            "license": str,
            "license_type": Or("permissive", "proprietary", "restricted"),
            "table_id": str,
            Optional(Regex(PROJECT_NAME_REGEX)): str,
            Optional("_frequency"): bool,
            Optional("_miranda_version"): bool,
            Optional("_remove_attrs"): Or(
                Regex(PROJECT_NAME_REGEX),
                {Regex(PROJECT_NAME_REGEX): Or(str, [str])},
            ),
            Optional(Regex(r"^_")): Or({str: str}, {str: {str: str}}),
        },
        _institution_in_header,
    )
}


variables_schema = {
    Optional("variables"): {
        Optional(str): {
            "standard_name": Regex(STANDARD_NAME_REGEX),
            "_cf_variable_name": str,
            Optional(*VARIABLE_TREATMENTS): Or(
                str, bool, {str: Or(str, bool, {str: str})}
            ),
            Optional("cell_methods"): Regex(CELL_METHODS_REGEX),
            Optional("comment"): str,
            Optional("description"): str,
            Optional("long_name"): str,
            Optional("original_long_name"): str,
            Optional("units"): str,
        }
    }
}


# Convert Schema
_convert_schema = Schema(
    And(
        dict,
        {
            **header_schema,
            **cf_dimensions_schema,
            **variables_schema,
        },
    ),
    ignore_extra_keys=False,  # Extra entries will raise a ValidationError
)


# This function accepts a path to a JSON file, loads the JSON data, and validates it using schema
def validate_json(json_file: str | Path, schema: Schema | None = None) -> bool:
    """
    Validate a JSON file against a schema.

    Parameters
    ----------
    json_file: str or pathlib.Path
        The path to the JSON file.
    schema: Schema, optional
        The schema to validate against.
        If None, will choose a definition based on filename parameters.

    Returns
    -------
    bool
        True if the JSON file is valid, False otherwise.
    """
    if not Path(json_file).is_file():
        msg = f"{json_file} is not a file."
        raise ValueError(msg)

    if schema is None:
        if "_cf_" in Path(json_file).name:
            schema = _convert_schema
        else:
            raise ValueError("No schema specified.")
    elif not isinstance(schema, Schema):
        raise ValueError("'schema' must be a Schema instance.")

    try:
        with Path(json_file).open() as f:
            data = json.load(f)
        _schema.validate(data)
        return True
    except (OSError, json.JSONDecodeError, SchemaError) as e:
        msg = f"Error validating JSON file {json_file}: {e}"
        logging.error(msg)
        return False
