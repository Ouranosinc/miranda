from __future__ import annotations
import logging
import warnings
from typing import Any

import numpy as np
import xarray as xr
from xclim.core.calendar import parse_offset

from miranda.treatments.utils import _get_section_entry_key, _iter_entry_key  # noqa
from miranda.units import check_time_frequency


__all__ = [
    "dimensions_compliance",
    "ensure_correct_time_frequency",
    "find_project_variable_codes",
    "offset_time_dimension",
]


def find_project_variable_codes(code: str, configuration: dict[str, Any]) -> str:
    """
    Find the variable code for a given variable name and project.

    Parameters
    ----------
    code : str
        Variable name.
    configuration : dict
        Configuration dictionary.

    Returns
    -------
    str
    """
    variable_codes = {}

    if "variables" not in configuration:
        raise ValueError("No `variables` section found in configuration. Check JSON.")

    for variable_code in configuration["variables"]:
        variable_name = configuration["variables"][variable_code].get("_variable_name")
        if variable_name:
            variable_codes[variable_name] = variable_code
        else:
            warnings.warn(
                f"Variable `{variable_code}` does not have accompanying `variable_name`. "
                f"Verify JSON. Continuing with `{variable_code}` as `variable_name`.",
                stacklevel=2,
            )
            variable_codes[variable_code] = variable_code

    if code in variable_codes.values():
        variable = code
    else:
        variable = variable_codes.get(code)
    if not variable:
        raise NotImplementedError(f"Variable `{code}` not supported.")

    return variable


def dimensions_compliance(ds: xr.Dataset, project: str, metadata: dict) -> xr.Dataset:
    """
    Rename dimensions to CF to their equivalents and reorder them if needed.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with dimensions to be updated.
    project : str
        Dataset project name.
    metadata : dict
        Metadata definition dictionary for project and variable(s).

    Returns
    -------
    xarray.Dataset
    """
    rename_dims = dict()
    for dim in ds.dims:
        if dim in metadata["dimensions"].keys():
            cf_name = _get_section_entry_key(metadata, "dimensions", dim, "_cf_dimension_name", project)
            if cf_name:
                rename_dims[dim] = cf_name

    # Rename dimensions
    _rename_dims = [str(d) for d in rename_dims.keys()]
    msg = f"Renaming dimensions: {', '.join(_rename_dims)}."
    logging.info(msg)
    ds = ds.rename(rename_dims)
    for new in ["lon", "lat"]:
        if new == "lon" and "lon" in ds.coords:
            if np.any(ds.lon > 180):
                lon1 = ds.lon.where(ds.lon <= 180.0, ds.lon - 360.0)
                ds[new] = lon1

        coord_precision = _get_section_entry_key(metadata, "dimensions", new, "_precision", project)
        if coord_precision is not None:
            ds[new] = ds[new].round(coord_precision)

    # Ensure that lon and lat are written in proper order for plotting purposes
    logging.info("Reordering dimensions.")
    transpose_order = []
    if "lat" in ds.dims and "lon" in ds.dims:
        transpose_order = ["lat", "lon"]
    elif "rlat" in ds.dims and "rlon" in ds.dims:
        transpose_order = ["rlat", "rlon"]
    if "time" in ds.dims and transpose_order:
        transpose_order.insert(0, "time")
        transpose_order.extend(list(set(ds.dims) - set(transpose_order)))
    ds = ds.transpose(*transpose_order)
    ds = ds.sortby(transpose_order)

    # Add dimension original name and update attrs
    logging.info("Updating dimension attributes.")
    dim_descriptions = metadata["dimensions"]
    for dim in metadata["dimensions"].keys():
        cf_name = dim_descriptions[dim].get("_cf_dimension_name")
        if cf_name is not None and cf_name in ds.dims:
            ds[cf_name].attrs.update(dict(original_variable=dim))
        else:
            # variable name already follows CF standards
            cf_name = dim
        for field in dim_descriptions[dim].keys():
            if not field.startswith("_"):
                ds[cf_name].attrs.update({field: dim_descriptions[dim][field]})

    prev_history = ds.attrs.get("history", "")
    history = f"Transposed and renamed dimensions. {prev_history}"
    ds.attrs.update(dict(history=history))

    return ds


def ensure_correct_time_frequency(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Ensure that time frequency is consistent with expected frequency for project."""
    key = "_ensure_correct_time"
    strict_time = "_strict_time"

    if "time" not in m["dimensions"].keys():
        msg = f"No time corrections listed for project `{p}`. Continuing..."
        warnings.warn(msg, stacklevel=2)
        return d

    if "time" not in list(d.variables.keys()):
        msg = f"No time dimension among data variables: {' ,'.join([str(v) for v in d.variables.keys()])}. Continuing..."
        logging.info(msg)
        return d

    if key in m["dimensions"]["time"].keys():
        freq_found = xr.infer_freq(d.time)
        if strict_time in m["dimensions"]["time"].keys():
            if not freq_found:
                msg = "Time frequency could not be found. There may be missing timesteps."
                if m["dimensions"]["time"].get(strict_time):
                    raise ValueError(msg)
                else:
                    warnings.warn(f"{msg} Continuing...", stacklevel=2)
                    return d

        correct_time_entry = m["dimensions"]["time"][key]
        if isinstance(correct_time_entry, str):
            correct_times = [parse_offset(correct_time_entry)[1]]
        elif isinstance(correct_time_entry, dict):
            correct_times = correct_time_entry.get(p)
            if isinstance(correct_times, list):
                correct_times = [parse_offset(t)[1] for t in correct_times]
            if correct_times is None:
                warnings.warn(f"No expected times set for specified project `{p}`.", stacklevel=2)
        elif isinstance(correct_time_entry, list):
            correct_times = correct_time_entry
        else:
            warnings.warn("No expected times set for family of projects.", stacklevel=2)
            return d

        if freq_found not in correct_times:
            error_msg = (
                f"Time frequency {freq_found} not among allowed frequencies: "
                f"{', '.join(correct_times) if isinstance(correct_times, list) else correct_times}"
            )
            if isinstance(correct_time_entry, dict):
                error_msg = f"{error_msg} for project `{p}`."
            else:
                error_msg = f"{error_msg}."
            raise ValueError(error_msg)

        msg = f"Resampling dataset with time frequency: {freq_found}."
        logging.info(msg)
        with xr.set_options(keep_attrs=True):
            d_out = d.assign_coords(time=d.time.resample(time=freq_found).mean(dim="time").time)
            d_out.time.attrs.update(d.time.attrs)

        prev_history = d.attrs.get("history", "")
        history = f"Resampled time with `freq={freq_found}`. {prev_history}"
        d_out.attrs.update(dict(history=history))
        return d_out

    return d


def offset_time_dimension(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Offset time dimension using listed frequency."""
    key = "_offset_time"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    converted = []
    offset, offset_meaning = None, None

    time_freq = dict()
    expected_period = _get_section_entry_key(m, "dimensions", "time", "_ensure_correct_time", p)
    if isinstance(expected_period, str):
        time_freq["expected_period"] = expected_period

    for vv, offs in _iter_entry_key(d, m, "dimensions", key, p):
        if offs:
            # Offset time by value of one time-step
            if offset is None and offset_meaning is None:
                try:
                    offset, offset_meaning = check_time_frequency(d, **time_freq)
                except TypeError:
                    msg = "Unable to parse the time frequency. Verify data integrity before retrying."
                    logging.error(msg)
                    raise

            msg = f"Offsetting data for `{vv}` by `{offset[0]} {offset_meaning}(s)`."
            logging.info(msg)
            with xr.set_options(keep_attrs=True):
                out = d[vv]
                out["time"] = out.time - np.timedelta64(offset[0], offset[1])
                d_out[vv] = out
            converted.append(vv)
        elif offs is False:
            msg = f"No time offsetting needed for `{vv}` in `{p}` (Explicitly set to False)."
            logging.info(msg)
            continue
        prev_history = d.attrs.get("history", "")
        history = f"Offset variable `{vv}` values by `{offset[0]} {offset_meaning}(s). {prev_history}"
        d_out.attrs.update(dict(history=history))

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]
    return d_out
