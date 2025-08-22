from __future__ import annotations

from schema import Optional, Or, Regex, Schema

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

cf_variables_schema = Schema(
    {
        str: {
            "standard_name": Regex(STANDARD_NAME_REGEX),
            "_cf_variable_name": str,
            Optional(Or(*VARIABLE_TREATMENTS)): Or(
                str, bool, Schema({str: Or(str, bool, Schema({str: str}))})
            ),
            Optional(Or(*CONTEXT_VARIABLE_TREATMENTS)): Or(
                Schema({str: Schema({str: Or(str, bool)})})
            ),
            Optional("cell_methods"): Regex(CELL_METHODS_REGEX),
            Optional("comments"): str,
            Optional("description"): str,
            Optional("long_name"): str,
            Optional("original_long_name"): str,
            Optional("units"): str,
        }
    },
    name="variables_schema",
)
