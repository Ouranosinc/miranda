from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any

_config_folder = Path(__file__).resolve().parent / "configs"


__all__ = ["load_json_data_mappings", "find_project_variable_codes"]


def load_json_data_mappings(project: str) -> dict[str, Any]:
    """Load JSON mappings for supported dataset conversions.

    Parameters
    ----------
    project : str

    Returns
    -------
    dict[str, Any]
    """
    if project == "eccc-homogenized":
        metadata_definition = json.load(
            open(_config_folder / "eccc-homogenized_attrs.json")
        )
    elif project == "eccc-obs":
        metadata_definition = json.load(open(_config_folder / "eccc-obs_attrs.json"))
    elif project == "eccc-obs-summary":
        metadata_definition = json.load(
            open(_config_folder / "eccc-obs-summary_attrs.json")
        )
    else:
        raise NotImplementedError(f"Project not supported: {project}")

    return metadata_definition


def find_project_variable_codes(code: str, table: str) -> str:
    """Find the variable code for a given variable name and project.

    Parameters
    ----------
    code : str
        Variable name.
    table : str
        Project name.

    Returns
    -------
    str
    """
    config = load_json_data_mappings(table)
    variable_codes = {}
    for variable_code in config["variables"]:
        variable_name = config["variables"][variable_code].get("_variable_name")
        if variable_name:
            variable_codes[variable_name] = variable_code
        else:
            warnings.warn(
                f"Variable `{variable_code}` does not have accompanying `variable_name`. "
                f"Verify JSON. Continuing with `{variable_code}` as `variable_name`."
            )
            variable_codes[variable_code] = variable_code

    if code in variable_codes.values():
        variable = code
    else:
        variable = variable_codes.get(code)
    if not variable:
        raise NotImplementedError(f"Variable `{code}` not supported.")

    return variable
