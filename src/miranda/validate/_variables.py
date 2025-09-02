from __future__ import annotations

from schema import And, Optional, Or, Regex, Schema

from ._regex import CELL_METHODS_REGEX, STANDARD_NAME_REGEX


__all__ = ["cf_variables_schema"]

VARIABLE_TREATMENTS = [
    "_corrected_standard_name",
    "_corrected_units",
    "_cf_units_conversion",
    "_correct_unit_names",
    "_invert_sign",
    "_offset_time",
    "_transformation",
    "_units_context",
    "_variable_conversion",
]

CONTEXT_VARIABLE_TREATMENTS = [
    "_clip_values",
]


def _standard_name_in_variable(variable_dict: dict):
    """
    Check for standard_name in the variable.

    Parameters
    ----------
    variable_dict : dict
        The Schema dictionary for the variable.

    Returns
    -------
    Schema
        The validated variable Schema dictionary.

    Raises
    ------
    ValueError
        If the variable dict has neither "standard_name" nor {"_standard_name": project_name} or {"_standard_name": false}
        If the variable dict contains both "standard_name" and {"_standard_name": project_name}
    """
    if not variable_dict or not isinstance(variable_dict, dict):
        raise ValueError("'time' must be present and be a dictionary.")

    # Must contain either 'source' or '_source', but not both
    has_standard_name = "standard_name" in variable_dict
    has_dynamic_standard_name = "_standard_name" in variable_dict

    if has_standard_name and has_dynamic_standard_name:  # both true or both false
        raise ValueError("Time dimension may contain either 'units' or '_units', but not both")

    if has_standard_name:
        if not isinstance(variable_dict["standard_name"], str):
            raise ValueError("'standard_name' must be a string")

    if has_dynamic_standard_name:
        dynamic_source = variable_dict["_standard_name"]
        if isinstance(dynamic_source, bool):
            if dynamic_source is True:
                raise ValueError("'_standard_name' cannot be True, must be a dict of {str: str} or False")

        elif not isinstance(dynamic_source, dict) or not any(isinstance(k, str) and isinstance(v, str) for k, v in dynamic_source.items()):
            raise ValueError("'_source' must be a dict of {str: str}")

    return variable_dict


cf_variables_schema = Schema(
    And(
        Schema(
            {
                str: {
                    Optional("standard_name"): Regex(STANDARD_NAME_REGEX),
                    "_cf_variable_name": str,
                    Optional(Or(*VARIABLE_TREATMENTS)): Or(str, bool, Schema({str: Or(str, bool, Schema({str: str}))})),
                    Optional(Or(*CONTEXT_VARIABLE_TREATMENTS)): Or(Schema({str: Schema({str: Or(str, bool)})})),
                    Optional("cell_methods"): Regex(CELL_METHODS_REGEX),
                    Optional("comments"): str,
                    Optional("description"): str,
                    Optional("long_name"): str,
                    Optional("original_long_name"): str,
                    Optional("units"): str,
                    Optional("_standard_name"): False,
                }
            },
            name="variables_schema",
        ),
        _standard_name_in_variable,
    )
)
