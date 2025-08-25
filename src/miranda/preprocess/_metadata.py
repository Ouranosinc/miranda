from __future__ import annotations
import logging
from typing import Any

from miranda import __version__ as __miranda_version__
from miranda.treatments.utils import load_json_data_mappings


__all__ = [
    "eccc_variable_metadata",
    "homogenized_column_definitions",
    "obs_column_definitions",
]


def eccc_variable_metadata(
    variable_code: str | int,
    project: str,
    generation: int | None = None,
    metadata: dict | None = None,
) -> dict[str, Any]:
    """
    Return the metadata for a given variable code and project.

    Parameters
    ----------
    variable_code: str or int
    project: {"eccc-ahccd", "eccc-obs", "eccc-obs-summary"}
    generation: {1, 2, 3}, optional
    metadata: dict, optional

    Returns
    -------
    dict
    """
    if project == "eccc-ahccd":
        generation = {1: "First", 2: "Second", 3: "Third"}.get(generation)
        if not generation:
            raise NotImplementedError(f"Generation '{generation}' not supported")
    else:
        generation = None

    if not metadata:
        metadata = load_json_data_mappings(project)

    if isinstance(variable_code, int):
        variable_code = str(variable_code).zfill(3)

    # code = find_project_variable_codes(variable_code, metadata)

    # Variable metadata
    variable_meta = metadata["variables"].get(variable_code)
    if variable_meta is None:
        raise ValueError(f"No metadata found for variable code: {variable_code}")

    variable_name = ""
    variable_name_fields = ["_variable_name", "_cf_variable_name"]
    if set(variable_name_fields).issubset(variable_meta.keys()):
        for variable_field in variable_name_fields:
            variable_name = variable_meta.get(variable_field)
            if variable_name:
                variable_meta["original_variable_code"] = variable_code
                del variable_meta[variable_field]
                variable_meta = {variable_name: variable_meta}
    else:
        variable_meta = {variable_code: variable_meta}
    if not variable_name:
        variable_name = variable_code

    # Dataset metadata
    header = metadata.get("Header")
    # Static handling of version global attributes
    miranda_version = header.get("_miranda_version")
    if miranda_version:
        if isinstance(miranda_version, bool):
            header["miranda_version"] = __miranda_version__
        elif isinstance(miranda_version, dict):
            if project in miranda_version.keys():
                header["miranda_version"] = __miranda_version__
        else:
            msg = f"`_miranda_version` not properly configured for project `{project}`. Not appending."
            logging.warning(msg)
    if "_miranda_version" in header:
        del header["_miranda_version"]

    to_delete = []
    # Conditional handling of global attributes based on fields
    for field in [f for f in header if f.startswith("_")]:
        if isinstance(header[field], bool):
            if header[field] and field == "_variable":
                header[field[1:]] = variable_name
        elif isinstance(header[field], dict) and generation:
            attr_treatment = header[field]["generation"]
            if field in ["_citation_product"]:
                for attribute, value in attr_treatment.items():
                    if attribute == generation:
                        header[field[1:]] = value
        else:
            raise AttributeError(f"Attribute treatment configuration for field `{field}` is not properly configured. Verify JSON.")
        to_delete.append(field)

    for field in to_delete:
        del header[field]

    return dict(metadata=variable_meta, header=header)


def homogenized_column_definitions(
    variable_code: str,
) -> tuple[dict, list[tuple[int, int]], dict[str, type[str | int | float] | Any], int]:
    """
    Return the column names, widths, and data types for the AHCCD fixed-width format data.

    Parameters
    ----------
    variable_code : str

    Returns
    -------
    tuple[dict, list[tuple[int, int]], dict[str, type[str | int | float] | Any], int]
    """
    metadata = load_json_data_mappings("eccc-homogenized")

    variable = metadata["variables"][variable_code]["_variable_name"]
    if variable.startswith("tas"):
        column_dtypes = {
            "No": str,
            "StnId": str,
            "Station name": str,
            "Prov": str,
            "FromYear": int,
            "FromMonth": int,
            "ToYear": int,
            "ToMonth": int,
            "%Miss": float,
            "Lat(deg)": float,
            "Long(deg)": float,
            "Elev(m)": int,
            "Joined": str,
            "RCS": str,
        }
        column_spaces = [(0, 5), (5, 6), (6, 8), (8, 9)]
        ii = 9
        # 31 days in a month
        for _ in range(31):
            column_spaces.append((ii, ii + 7))
            ii += 7
            column_spaces.append((ii, ii + 1))
            ii += 1
        header_row = 3

    elif variable.startswith("pr"):
        column_dtypes = {
            "Prov": str,
            "Station name": str,
            "stnid": str,
            "beg yr": int,
            "beg mon": int,
            "end yr": int,
            "end mon": int,
            "lat (deg)": float,
            "long (deg)": float,
            "elev (m)": int,
            "stns joined": str,
        }
        column_spaces = [(0, 4), (4, 5), (5, 7), (7, 8)]
        ii = 8
        # 31 days in a month
        for _ in range(31):
            column_spaces.append((ii, ii + 8))
            ii += 8
            column_spaces.append((ii, ii + 1))
            ii += 1
        header_row = 0

    else:
        raise KeyError

    column_names = {col.lower().split("(")[0].replace("%", "pct_").strip().replace(" ", "_"): col for col in list(column_dtypes.keys())}

    return column_names, column_spaces, column_dtypes, header_row


def obs_column_definitions(
    time_frequency: str,
) -> tuple[list[str], list[int], list[type[str | int]], int]:
    """Return the column names, widths, and data types for the fixed-width format."""
    if time_frequency.lower() in ["h", "hour", "hourly"]:
        num_observations = 24
        column_names = ["code", "year", "month", "day", "code_var"]
        column_widths = [7, 4, 2, 2, 3]
        column_dtypes = [str, int, int, int, str]
    elif time_frequency.lower() in ["d", "day", "daily"]:
        num_observations = 31
        column_names = ["code", "year", "month", "code_var"]
        column_widths = [7, 4, 2, 3]
        column_dtypes = [str, int, int, str]
    else:
        raise NotImplementedError("`mode` must be 'h'/'hourly or 'd'/'daily'.")

    header = 0

    # Add the data columns
    for i in range(1, num_observations + 1):
        data_entry, flag_entry = f"D{i:0n}", f"F{i:0n}"
        column_names.append(data_entry)
        column_names.append(flag_entry)
        column_widths.extend([6, 1] * num_observations)
        column_dtypes.extend([str, str])

    return column_names, column_widths, column_dtypes, header
