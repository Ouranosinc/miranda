from __future__ import annotations

import json
import logging
from pathlib import Path

from schema import And, Optional, Or, Regex, Schema, SchemaError

__all__ = ["validate_json"]


PROJECT_NAME_REGEX = r"^[a-zA-Z]\S+"
VALID_TIME_FREQUENCIES = r"^(\d+)?(Y|YS|M|MS|W|D|H|T|min|S|L|ms|U|us|N)$"

# Example schema
_default = Schema(
    {
        "Header": {
            "Conventions": Regex(r"CF-\d\.\d+"),
            "source": str,
            "type": str,
            "institution": str,
            "processing_level": Or("raw", "biasadjusted"),
            "license": str,
            "license_type": Or("permissive", "proprietary"),
            "table_id": str,
            Optional(Regex(PROJECT_NAME_REGEX)): str,
            Optional("_frequency"): bool,
            Optional("_miranda_version"): bool,
            Optional("_remove_attrs"): Or(
                Regex(PROJECT_NAME_REGEX), {Regex(PROJECT_NAME_REGEX): Or(str, [str])}
            ),
            Optional(Regex(r"^_")): Or({str: str}, {str: {str: str}}),
        },
        Optional("dimensions"): {
            Optional("time"): {
                Optional("_ensure_correct_time"): Or(
                    Regex(VALID_TIME_FREQUENCIES),
                    {Regex(PROJECT_NAME_REGEX): Regex(VALID_TIME_FREQUENCIES)},
                )
            }
        },
        "variables": dict,
    }
)


# This function accepts a path to a JSON file, loads the JSON data, and validates it using schema
def validate_json(json_file: str | Path, _schema: Schema = _default) -> bool:
    """
    Validate a JSON file against a schema.

    Parameters
    ----------
    json_file: str or pathlib.Path
        The path to the JSON file.
    _schema: Schema
        The schema to validate against.
        Default is a simple schema that checks for name, age, and optional email and address.

    Returns
    -------
    bool
        True if the JSON file is valid, False otherwise.
    """
    try:
        with Path(json_file).open() as f:
            data = json.load(f)
        _schema.validate(data)
        return True
    except (OSError, json.JSONDecodeError, SchemaError) as e:
        msg = f"Error validating JSON file {json_file}: {e}"
        logging.error(msg)
        return False
