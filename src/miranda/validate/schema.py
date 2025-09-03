"""Validate outputted metadata against CF-like schemas."""

from __future__ import annotations
import json
import logging
from pathlib import Path

from schema import And, Optional, Or, Regex, Schema, SchemaError

from ._dimensions import cf_dimensions_schema
from ._regex import (
    CF_CONVENTIONS_REGEX,
    PROJECT_NAME_REGEX,
)
from ._variables import cf_variables_schema


__all__ = [
    "cf_dimensions_schema",
    "cf_header_schema",
    "cf_variables_schema",
    "converter_schema",
    "validate_json",
]

LICENSES_TYPES = ["open", "permissive", "proprietary", "restricted"]


def _source_in_header(header_dict: dict):
    """
    Check for source in the header

    Parameters
    ----------
    header_dict : dict
        The Schema dictionary for the header.

    Returns
    -------
    Schema
        The validated header Schema dictionary.

    Raises
    ------
    ValueError
        If the header dict has neither "source" nor {"_source": project_name}
        If the header dict contains both "source" and {"_source": project_name}
    """
    if not header_dict or not isinstance(header_dict, dict):
        raise ValueError("'time' must be present and be a dictionary.")

    # Must contain either 'source' or '_source', but not both
    has_source = "source" in header_dict
    has_dynamic_source = "_source" in header_dict

    if has_source and has_dynamic_source:  # both true or both false
        raise ValueError("Time dimension may contain either 'units' or '_units', but not both")

    if has_source:
        if not isinstance(header_dict["source"], str):
            raise ValueError("'source' must be a string")

    if has_dynamic_source:
        dynamic_source = header_dict["_source"]
        if not isinstance(dynamic_source, dict) or not any(isinstance(k, str) and isinstance(v, str) for k, v in dynamic_source.items()):
            raise ValueError("'_source' must be a dict of {str: str}")

    return header_dict


def _institution_in_header(header_dict: dict):
    """
    Check for institution metadata in header.

    Parameters
    ----------
    header_dict : dict
        The Schema dictionary for the header.

    Returns
    -------
    Schema
        The validated header Schema dictionary.

    Raises
    ------
    ValueError
        If the Header has neither "institution" nor {"_map_attrs": {str : "institution"}
        If the Header contains both "institution" and {"_map_attrs": {str : "institution"}
    """
    # Must contain either 'institution' or '_map_attrs', but not both
    has_institution = "institution" in header_dict
    has_map_attrs = "_map_attrs" in header_dict

    if has_institution == has_map_attrs:  # both true or both false
        raise ValueError("Header must contain either 'institution' or '_map_attrs', but not both")

    if has_institution:
        if not isinstance(header_dict["institution"], str):
            raise ValueError("'institution' must be a string")

    if has_map_attrs:
        map_attrs = header_dict["_map_attrs"]
        if not isinstance(map_attrs, dict) or not any(isinstance(k, str) and v == "institution" for k, v in map_attrs.items()):
            raise ValueError("'_map_attrs' must be a dict of {str: 'institution'}")

    return header_dict


cf_header_schema = Schema(
    And(
        Schema(
            {
                "Conventions": Regex(CF_CONVENTIONS_REGEX),
                "type": str,  # FIXME: Should this be constrained to specific values? e.g., "simulation", "observation", etc.
                "processing_level": Or("raw", "biasadjusted"),
                Optional(Regex(r"^license$|^licence$")): str,
                Regex(r"^license_type$|^licence_type$"): Or(
                    *LICENSES_TYPES,
                    Schema({Regex(PROJECT_NAME_REGEX): Or(*LICENSES_TYPES)}),
                ),
                "table_id": str,
                Optional(Regex(PROJECT_NAME_REGEX)): str,
                Optional("_frequency"): bool,
                Optional(Regex(r"^_license$|^_licence$")): {str: Or(str, Schema({str: str}))},
                Optional("_miranda_version"): bool,
                Optional("_remove_attrs"): Or(
                    Schema(Regex(PROJECT_NAME_REGEX)),
                    Schema({Regex(PROJECT_NAME_REGEX): Or(str, Schema([str]))}),
                ),
                Optional(Regex(r"^_")): {str: Or(str, bool, Schema({str: str}))},
            }
        ),
        _institution_in_header,
        _source_in_header,
    ),
    name="header_schema",
)


# Converter Schema
converter_schema = Schema(
    {
        "Header": cf_header_schema,
        "variables": cf_variables_schema,
        "dimensions": cf_dimensions_schema,
    },
    ignore_extra_keys=False,  # Extra entries will raise a ValidationError
    name="convert_schema",
)


# This function accepts a path to a JSON file, loads the JSON data, and validates it using schema
def validate_json(json_file: str | Path, schema: Schema | None = None) -> bool:
    """
    Validate a JSON file against a schema.

    Parameters
    ----------
    json_file : str or pathlib.Path
        The path to the JSON file.
    schema : Schema, optional
        The schema to validate against.
        If None, will choose a definition based on filename parameters.

    Returns
    -------
    bool
        True if the JSON file is valid, False otherwise.

    Raises
    ------
    ValueError
        If the JSON file does not exist, or if the schema is not CF-compliant.
    OSError
        If there is an error reading the JSON file.
    json.JSONDecodeError
        If the JSON file is not valid JSON.
    SchemaError
        If the JSON data does not conform to the schema.
    """
    if not Path(json_file).is_file():
        msg = f"{json_file} is not a file."
        raise ValueError(msg)

    if schema is None:
        if "_cf_" in Path(json_file).name:
            schema = converter_schema
        else:
            raise ValueError("Schema is not CF-compliant. No validation is possible.")
    elif not isinstance(schema, Schema):
        raise ValueError("'schema' must be a Schema instance.")

    try:
        with Path(json_file).open() as f:
            data = json.load(f)
        schema.validate(data)
        return True
    except (OSError, json.JSONDecodeError) as e:
        msg = f"Error validating JSON file {json_file}: {e}"
        logging.error(msg)
        raise
    except SchemaError as e:
        msg = f"Schema validation error in {json_file}: {e}"
        logging.error(msg)
        raise
