from __future__ import annotations

import json
import logging
from pathlib import Path

from schema import And, Optional, Regex, Schema, SchemaError, Or

__all__ = ["validate_json"]

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
            Optional(Regex("^[a-zA-Z]")): str,
            Optional("_frequency"): bool,
            Optional("_miranda_version"): bool,
            Optional(Regex(r"^_")): {str: str},

        },
        "dimensions": dict,
        "variables": dict
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
        with open(json_file) as f:
            data = json.load(f)
        _schema.validate(data)
        return True
    except (OSError, json.JSONDecodeError, SchemaError) as e:
        msg = f"Error validating JSON file {json_file}: {e}"
        logging.error(msg)
        return False
