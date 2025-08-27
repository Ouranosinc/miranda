"""Validation utilities and definitions."""

__all__ = [
    "CELL_METHODS_REGEX",
    "CF_CONVENTIONS_REGEX",
    "PROJECT_NAME_REGEX",
    "STANDARD_NAME_REGEX",
    "TIME_UNITS_REGEX",
    "VALID_TIME_FREQUENCY_REGEX",
]

CELL_METHODS_REGEX = r"^\w+:\s*\w+"
CF_CONVENTIONS_REGEX = r"CF-\d\.\d+"
PROJECT_NAME_REGEX = r"^[a-zA-Z]\S+"
STANDARD_NAME_REGEX = r"^\w+(?:_\w+)*$"
TIME_UNITS_REGEX = r"^(seconds|minutes|hours|days|months|years)\s+since\s+(\d{4})-(\d{1,2})-(\d{1,2})"

VALID_TIME_FREQUENCY_REGEX = r"^(\d+)?(Y|YS|Q|QS|M|MS|W|D|h|m|s|ms|us|ns)$"
