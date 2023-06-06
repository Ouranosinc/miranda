from __future__ import annotations

import datetime
import json
import logging.config
import os
import warnings
from collections.abc import Iterator, Sequence
from functools import partial
from pathlib import Path
from typing import Any, Callable

import numpy as np
import xarray as xr
import xclim.core.units
from xarray.coding import times
from xclim.core import units
from xclim.core.calendar import parse_offset

from miranda import __version__ as __miranda_version__
from miranda.gis import subset_domain
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

from .utils import date_parser, find_version_hash

logging.config.dictConfig(LOGGING_CONFIG)

VERSION = datetime.datetime.now().strftime("%Y.%m.%d")

__all__ = [
    "dataset_corrections",
    "dims_conversion",
    "dataset_conversion",
    "load_json_data_mappings",
    "metadata_conversion",
    "threshold_mask",
    "variable_conversion",
]


def load_json_data_mappings(project: str) -> dict[str, Any]:
    """Load JSON mappings for supported dataset conversions.

    Parameters
    ----------
    project : str

    Returns
    -------
    dict[str, Any]
    """
    data_folder = Path(__file__).resolve().parent / "data"

    if project.startswith("era5"):
        metadata_definition = json.load(open(data_folder / "ecmwf_cf_attrs.json"))
    elif project in ["rdrs-v21"]:
        metadata_definition = json.load(open(data_folder / "eccc_rdrs_cf_attrs.json"))
    elif project in ["agcfsr", "agmerra2"]:  # This should handle the AG versions:
        metadata_definition = json.load(open(data_folder / "nasa_cf_attrs.json"))
    elif project in ["cordex", "cmip5", "cmip6"]:
        metadata_definition = json.load(open(data_folder / "cmip_ouranos_attrs.json"))
    elif project == "ets-grnch":
        metadata_definition = json.load(open(data_folder / "ets_grnch_cf_attrs.json"))
    elif project == "nrcan-gridded-10km":
        raise NotImplementedError()
    elif project == "wfdei-gem-capa":
        metadata_definition = json.load(open(data_folder / "usask_cf_attrs.json"))
    elif project.startswith("melcc"):
        metadata_definition = json.load(open(data_folder / "melcc_cf_attrs.json"))
    elif project.startswith("ec"):
        metadata_definition = json.load(open(data_folder / "eccc_cf_attrs.json"))
    elif project in ["NEX-GDDP-CMIP6"]:
        metadata_definition = json.load(open(data_folder / "nex-gddp-cmip6_attrs.json"))
    elif project in ["ESPO-G6-R2"]:
        metadata_definition = json.load(open(data_folder / "espo-g6-r2_attrs.json"))
    elif project in ["ESPO-G6-E5L"]:
        metadata_definition = json.load(open(data_folder / "espo-g6-e5l_attrs.json"))
    elif project in ["EMDNA"]:
        metadata_definition = json.load(open(data_folder / "emdna_cf_attrs.json"))
    else:
        raise NotImplementedError()

    return metadata_definition


def _get_section_entry_key(meta, entry, var, key, project):
    var_meta = meta[entry].get(var, {})
    if key in var_meta:
        if isinstance(var_meta[key], dict):
            config = var_meta[key].get(project)
            if config is None and "all" in var_meta[key].keys():
                config = var_meta[key].get("all")
            return config
        return var_meta[key]
    return None


def _iter_entry_key(ds, meta, entry, key, project):
    for vv in set(ds.data_vars).intersection(meta[entry]):
        val = _get_section_entry_key(meta, entry, vv, key, project)
        yield vv, val


def _simple_fix_dims(d: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
    """Adjust dimensions found in a file so that it can be used for regridding purposes."""
    if "lon" not in d.dims or "lat" not in d.dims:
        dim_rename = dict()
        for dim in d.dims:
            if str(dim).lower().startswith("lon"):
                dim_rename[str(dim)] = "lon"
            if str(dim).lower().startswith("lat"):
                dim_rename[str(dim)] = "lat"
        d = d.rename(dim_rename)
    if np.any(d.lon > 180):
        lon_wrapped = d.lon.where(d.lon <= 180.0, d.lon - 360.0)
        d["lon"] = lon_wrapped
        d = d.sortby(["lon"])

    if "time" in d.dims:
        d = d.isel(time=0, drop=True)

    return d


def conservative_regrid(
    ds: xr.DataArray | xr.Dataset, ref_grid: xr.DataArray | xr.Dataset
) -> xr.DataArray | xr.Dataset:
    """Perform a conservative_normed regridding"""
    try:
        import xesmf as xe  # noqa
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "This function requires the `xesmf` library which is not installed. "
            "Regridding step will be skipped."
        )

    ref_grid = _simple_fix_dims(ref_grid)
    method = "conservative_normed"

    logging.info(
        f"Performing regridding and masking with `xesmf` using method: {method}."
    )

    regridder = xe.Regridder(ds, ref_grid, method, periodic=False)
    ds = regridder(ds)

    ds.attrs["history"] = (
        f"{datetime.datetime.now()}:"
        f"Regridded dataset using xesmf with method: {method}. "
        f"{ds.attrs.get('history')}".strip()
    )
    return ds


def threshold_mask(
    ds: xr.Dataset | xr.DataArray,
    *,
    mask: xr.Dataset | xr.DataArray,
    mask_cutoff: float | bool = False,
) -> xr.Dataset | xr.DataArray:
    """Land-Sea mask operations.

    Parameters
    ----------
    ds : xr.Dataset or str or os.PathLike
    mask : xr.Dataset or xr.DataArray
    mask_cutoff : float or bool

    Returns
    -------
    xr.Dataset or xr.DataArray
    """
    mask = _simple_fix_dims(mask)

    if isinstance(mask, xr.Dataset):
        if len(mask.data_vars) == 1:
            mask_variable = list(mask.data_vars)[0]
            mask = mask[mask_variable]
        else:
            raise ValueError(
                "More than one data variable found in land-sea mask. Supply a DataArray instead."
            )
    else:
        mask_variable = mask.name

    try:
        from clisops.core import subset_bbox  # noqa

        log_msg = f"Masking dataset with {mask_variable}."
        if mask_cutoff:
            log_msg = f"{log_msg.strip('.')} at `{mask_cutoff}` cutoff value."
        logging.info(log_msg)

        lon_bounds = np.array([ds.lon.min(), ds.lon.max()])
        lat_bounds = np.array([ds.lat.min(), ds.lat.max()])

        mask_subset = subset_bbox(
            mask,
            lon_bnds=lon_bounds,
            lat_bnds=lat_bounds,
        ).load()
    except ModuleNotFoundError:
        log_msg = (
            "This function requires the `clisops` library which is not installed. "
            "subsetting step will be skipped."
        )
        warnings.warn(log_msg)
        mask_subset = mask.load()

    if mask_subset.dtype == bool:
        if mask_cutoff:
            logging.warning("Mask value cutoff set for boolean mask. Ignoring.")
        mask_subset = mask_subset.where(mask)
    else:
        mask_subset = mask_subset.where(mask >= mask_cutoff)
    ds = ds.where(mask_subset.notnull())

    if mask_subset.min() >= 0:
        if mask_subset.max() <= 1.00000001:
            cutoff_info = f"{mask_cutoff * 100} %"
        elif mask_subset.max() <= 100.00000001:
            cutoff_info = f"{mask_cutoff} %"
        else:
            cutoff_info = f"{mask_cutoff}"
    else:
        cutoff_info = f"{mask_cutoff}"
    ds.attrs["mask_cutoff"] = cutoff_info

    prev_history = ds.attrs.get("history", "")
    history_msg = f"Mask calculated using `{mask_variable}`."
    if mask_cutoff:
        history_msg = f"{history_msg.strip('.')} with cutoff value `{cutoff_info}`."
    history = f"{history_msg} {prev_history}".strip()
    ds.attrs.update(dict(history=history))

    return ds


def correct_time_entries(
    d: xr.Dataset,
    split: str = "_",
    location: int = -1,
    field: str = "time",
) -> xr.Dataset:
    filename = d.encoding["source"]
    date = date_parser(Path(filename).stem.split(split)[location])
    vals = np.arange(len(d[field]))
    days_since = f"days since {date}"
    time = xr.coding.times.decode_cf_datetime(
        vals, units=days_since, calendar="standard"
    )
    d = d.assign_coords({field: time})

    prev_history = d.attrs.get("history", "")
    history = (
        f"Time index recalculated in preprocessing step ({days_since}). {prev_history}"
    )
    d.attrs.update(dict(history=history))

    return d


def correct_var_names(d: xr.Dataset, split: str = "_", location: int = 0) -> xr.Dataset:
    """

    Parameters
    ----------
    d : xarray.Dataset
    split : str
    location : int

    Returns
    -------
    xarray.Dataset
    """
    filename = d.encoding["source"]
    new_name = Path(filename).stem.split(split)[location]
    old_name = list(d.data_vars.keys())[0]

    prev_history = d.attrs.get("history", "")
    history = f"Variable renamed in preprocessing step ({old_name}: {new_name}). {prev_history}"
    d.attrs.update(dict(history=history))

    return d.rename({old_name: new_name})


def preprocessing_corrections(ds: xr.Dataset, project: str) -> xr.Dataset:
    """Corrections function dispatcher to ensure minimal dataset validity on open.

    Parameters
    ----------
    ds : xarray.Dataset
    project : str

    Returns
    -------
    xarray.Dataset
    """

    def _preprocess_correct(d: xr.Dataset, *, ops: list[partial]) -> xr.Dataset:
        for correction in ops:
            d = correction(d)
        return d

    correction_fields = load_json_data_mappings(project).get("_preprocess")
    if correction_fields:
        preprocess_ops = []
        for field in correction_fields:
            if field == "_variable_name":
                preprocess_ops.append(
                    partial(correct_var_names, **correction_fields[field])
                )
            if field == "_time":
                preprocess_ops.append(
                    partial(correct_time_entries, **correction_fields[field])
                )
        if preprocess_ops:
            corrector = partial(_preprocess_correct, ops=preprocess_ops)
            return corrector(ds)
    return ds


def _correct_units_names(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    key = "_corrected_units"
    for var, val in _iter_entry_key(d, m, "variables", key, p):
        if val:
            d[var].attrs["units"] = val

    # FIXME: This is no longer relevant. Performed under dimension conversion step.
    val_time = _get_section_entry_key(m, "variables", "time", key, p)
    if val_time:
        d["time"].attrs["units"] = val_time

    return d


# for de-accumulation or conversion to flux
def _transform(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    key = "_transformation"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    converted = []
    offset, offset_meaning = None, None

    time_freq = dict()
    expected_period = _get_section_entry_key(
        m, "dimensions", "time", "_ensure_correct_time", p
    )
    if isinstance(expected_period, str):
        time_freq["expected_period"] = expected_period

    for vv, trans in _iter_entry_key(d, m, "variables", key, p):
        if trans:
            if trans == "deaccumulate":
                # Time-step accumulated total to time-based flux (de-accumulation)
                if offset is None and offset_meaning is None:
                    try:
                        offset, offset_meaning = get_time_frequency(d, **time_freq)
                    except TypeError:
                        logging.error(
                            "Unable to parse the time frequency. Verify data integrity before retrying."
                        )
                        raise

                logging.info(f"De-accumulating units for variable `{vv}`.")
                with xr.set_options(keep_attrs=True):
                    out = d[vv].diff(dim="time")
                    out = d[vv].where(
                        getattr(d[vv].time.dt, offset_meaning) == offset[0],
                        out.broadcast_like(d[vv]),
                    )
                    out = units.amount2rate(out, out_units=m["variables"][vv]["units"])
                    d_out[vv] = out
                converted.append(vv)
            elif trans == "amount2rate":
                # NOTE: This treatment is no longer needed in xclim v0.43.0+ but is kept for backwards compatibility
                # frequency-based totals to time-based flux
                logging.info(
                    f"Performing amount-to-rate units conversion for variable `{vv}`."
                )
                with xr.set_options(keep_attrs=True):
                    out = units.amount2rate(
                        d[vv],
                        out_units=m["variables"][vv]["units"],
                    )
                    d_out[vv] = out
                converted.append(vv)
            elif isinstance(trans, str):
                if trans.startswith("op "):
                    op = trans[3]
                    value = trans[4:].strip()
                    if value.startswith("attrs"):
                        value = units.str2pint(d[vv].attrs[value[6:]])
                    else:
                        value = units.str2pint(value)
                    with xr.set_options(keep_attrs=True):
                        if op == "+":
                            value = units.convert_units_to(value, d[vv])
                            d_out[vv] = d[vv] + value
                        elif op == "-":
                            value = units.convert_units_to(value, d[vv])
                            d_out[vv] = d[vv] - value
                        elif op == "*":
                            d_out[vv] = units.pint_multiply(d[vv], value)
                        elif op == "/":
                            d_out[vv] = units.pint_multiply(d[vv], 1 / value)
                        else:
                            raise NotImplementedError(
                                f"Op transform doesn't implement the «{op}» operator."
                            )
                converted.append(vv)
            else:
                raise NotImplementedError(f"Unknown transformation: {trans}")
        elif trans is False:
            logging.info(
                f"No transformations needed for `{vv}` (Explicitly set to False)."
            )
            continue

        prev_history = d.attrs.get("history", "")
        history = (
            f"Transformed variable `{vv}` values using method `{trans}`. {prev_history}"
        )
        d_out.attrs.update(dict(history=history))

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]
    return d_out


def _offset_time(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    key = "_offset_time"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    converted = []
    offset, offset_meaning = None, None

    time_freq = dict()
    expected_period = _get_section_entry_key(
        m, "dimensions", "time", "_ensure_correct_time", p
    )
    if isinstance(expected_period, str):
        time_freq["expected_period"] = expected_period

    for vv, offs in _iter_entry_key(d, m, "dimensions", key, p):
        if offs:
            # Offset time by value of one time-step
            if offset is None and offset_meaning is None:
                try:
                    offset, offset_meaning = get_time_frequency(d, **time_freq)
                except TypeError:
                    logging.error(
                        "Unable to parse the time frequency. Verify data integrity before retrying."
                    )
                    raise

            logging.info(
                f"Offsetting data for `{vv}` by `{offset[0]} {offset_meaning}(s)`."
            )
            with xr.set_options(keep_attrs=True):
                out = d[vv]
                out["time"] = out.time - np.timedelta64(offset[0], offset[1])
                d_out[vv] = out
            converted.append(vv)
        elif offs is False:
            logging.info(
                f"No time offsetting needed for `{vv}` in `{p}` (Explicitly set to False)."
            )
            continue
        prev_history = d.attrs.get("history", "")
        history = f"Offset variable `{vv}` values by `{offset[0]} {offset_meaning}(s). {prev_history}"
        d_out.attrs.update(dict(history=history))

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]
    return d_out


def _invert_sign(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    key = "_invert_sign"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    converted = []
    for vv, inv_sign in _iter_entry_key(d, m, "variables", key, p):
        if inv_sign:
            logging.info(f"Inverting sign for `{vv}` (switching direction of values).")
            with xr.set_options(keep_attrs=True):
                out = d[vv]
                d_out[out.name] = -out
            converted.append(vv)
        elif inv_sign is False:
            logging.info(
                f"No sign inversion needed for `{vv}` in `{p}` (Explicitly set to False)."
            )
            continue
        prev_history = d.attrs.get("history", "")
        history = f"Inverted sign for variable `{vv}` (switched direction of values). {prev_history}"
        d_out.attrs.update(dict(history=history))

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]
    return d_out


# For converting variable units to standard workflow units
def _units_cf_conversion(d: xr.Dataset, m: dict) -> xr.Dataset:
    if "time" in m["dimensions"].keys():
        if m["dimensions"]["time"].get("units"):
            d["time"]["units"] = m["dimensions"]["time"]["units"]

    for vv, unit in _iter_entry_key(d, m, "variables", "units", None):
        if unit:
            with xr.set_options(keep_attrs=True):
                d[vv] = units.convert_units_to(d[vv], unit, context="hydro")
            prev_history = d.attrs.get("history", "")
            history = f"Converted variable `{vv}` to CF-compliant units (`{unit}`). {prev_history}"
            d.attrs.update(dict(history=history))

    return d


# For clipping variable values to an established maximum/minimum
def _clip_values(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    key = "_clip_values"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    converted = []
    for vv in d.data_vars:
        if vv in m["variables"].keys():
            clip_values = _get_section_entry_key(m, "variables", vv, key, p)
            if clip_values:
                min_value, max_value = None, None
                # Gather unit conversion context, if applicable
                context = clip_values.get("context", None)
                for op, value in clip_values.items():
                    if op == "min":
                        min_value = xclim.core.units.convert_units_to(
                            value, d[vv], context
                        )
                    if op == "max":
                        max_value = xclim.core.units.convert_units_to(
                            value, d[vv], context
                        )
                logging.info(
                    f"Clipping min/max values for `{vv}` ({min_value}/{max_value})."
                )
                with xr.set_options(keep_attrs=True):
                    out = d[vv]
                    d_out[out.name] = out.clip(min_value, max_value)
                converted.append(vv)
            elif clip_values is False:
                logging.info(
                    f"No clipping of values needed for `{vv}` in `{p}` (Explicitly set to False)."
                )
                continue
            else:
                logging.info(f"No clipping of values needed for `{vv}` in `{p}`.")
                continue

            prev_history = d.attrs.get("history", "")
            history = f"Clipped variable `{vv}` with `min={min_value}` and `max={max_value}`. {prev_history}"
            d_out.attrs.update(dict(history=history))

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]

    return d_out


def _ensure_correct_time(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    key = "_ensure_correct_time"
    strict_time = "_strict_time"

    if "time" not in m["dimensions"].keys():
        logging.warning(f"No time corrections listed for project `{p}`. Continuing...")
        return d

    if "time" not in list(d.variables.keys()):
        logging.info(
            "No time dimension among data variables: "
            f"{' ,'.join([str(v) for v in d.variables.keys()])}. "
            "Continuing..."
        )
        return d

    if key in m["dimensions"]["time"].keys():
        freq_found = xr.infer_freq(d.time)
        if strict_time in m["dimensions"]["time"].keys():
            if not freq_found:
                msg = (
                    "Time frequency could not be found. There may be missing timesteps."
                )
                if m["dimensions"]["time"].get(strict_time):
                    raise ValueError(msg)
                else:
                    logging.warning(f"{msg} Continuing...")
                    return d

        correct_time_entry = m["dimensions"]["time"][key]
        if isinstance(correct_time_entry, str):
            correct_times = [parse_offset(correct_time_entry)[1]]
        elif isinstance(correct_time_entry, dict):
            correct_times = correct_time_entry.get(p)
            if isinstance(correct_times, list):
                correct_times = [parse_offset(t)[1] for t in correct_times]
            if correct_times is None:
                logging.warning(f"No expected times set for specified project `{p}`.")
        elif isinstance(correct_time_entry, list):
            correct_times = correct_time_entry
        else:
            logging.warning("No expected times set for family of projects.")
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

        logging.info(f"Resampling dataset with time frequency: {freq_found}.")
        with xr.set_options(keep_attrs=True):
            d_out = d.assign_coords(
                time=d.time.resample(time=freq_found).mean(dim="time").time
            )
            d_out.time.attrs.update(d.time.attrs)

        prev_history = d.attrs.get("history", "")
        history = f"Resampled time with `freq={freq_found}`. {prev_history}"
        d_out.attrs.update(dict(history=history))
        return d_out

    return d


# For renaming and reordering lat and lon dims
def dims_conversion(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Rename dimensions to CF to their equivalents.

    Parameters
    ----------
    d : xarray.Dataset
        Dataset with dimensions to be updated.
    p : str
        Dataset project name.
    m : dict
        Metadata definition dictionary for project and variable(s).

    Returns
    -------
    xarray.Dataset
    """
    rename_dims = dict()
    for dim in d.dims:
        if dim in m["dimensions"].keys():
            cf_name = _get_section_entry_key(
                m, "dimensions", dim, "_cf_dimension_name", p
            )
            if cf_name:
                rename_dims[dim] = cf_name
    d = d.rename(rename_dims)
    for new in ["lon", "lat"]:
        if new == "lon" and "lon" in d.coords:
            if np.any(d.lon > 180):
                lon1 = d.lon.where(d.lon <= 180.0, d.lon - 360.0)
                d[new] = lon1

        coord_precision = _get_section_entry_key(m, "dimensions", new, "_precision", p)
        if coord_precision is not None:
            d[new] = d[new].round(coord_precision)

    # Ensure that lon and lat are written in proper order for plotting purposes
    transpose_order = []
    if "lat" in d.dims and "lon" in d.dims:
        transpose_order = ["lat", "lon"]
    elif "rlat" in d.dims and "rlon" in d.dims:
        transpose_order = ["rlat", "rlon"]
    if "time" in d.dims and transpose_order:
        transpose_order.insert(0, "time")
        transpose_order.extend(list(set(d.dims) - set(transpose_order)))
    d = d.transpose(*transpose_order)
    d = d.sortby(transpose_order)

    # Add dimension original name and update attrs
    dim_descriptions = m["dimensions"]
    for dim in m["dimensions"].keys():
        cf_name = dim_descriptions[dim].get("_cf_dimension_name")
        if cf_name is not None and cf_name in d.dims:
            d[cf_name].attrs.update(dict(original_variable=dim))
        else:
            # variable name already follows CF standards
            cf_name = dim
        for field in dim_descriptions[dim].keys():
            if not field.startswith("_"):
                d[cf_name].attrs.update({field: dim_descriptions[dim][field]})

    prev_history = d.attrs.get("history", "")
    history = f"Transposed and renamed dimensions. {prev_history}"
    d.attrs.update(dict(history=history))

    return d


def variable_conversion(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Add variable metadata and remove nonstandard entries.

    Parameters
    ----------
    d : xarray.Dataset
        Dataset with variable(s) to be updated.
    p : str
        Dataset project name.
    m : dict
        Metadata definition dictionary for project and variable(s).

    Returns
    -------
    xarray.Dataset
    """
    var_descriptions = m["variables"]
    var_correction_fields = [
        "_clip_values",
        "_corrected_units",
        "_invert_sign",
        "_offset_time",
        "_transformation",
    ]
    for var in d.variables:
        if var in var_descriptions.keys():
            for field in var_correction_fields:
                if field in var_descriptions[var].keys():
                    del var_descriptions[var][field]
            d[var].attrs.update(var_descriptions[var])

    # Rename data variables
    for orig_var_name, cf_name in _iter_entry_key(
        d, m, "variables", "_cf_variable_name", None
    ):
        if cf_name is not None:
            d = d.rename({orig_var_name: cf_name})
            d[cf_name].attrs.update(dict(original_variable=orig_var_name))
            del d[cf_name].attrs["_cf_variable_name"]

    return d


def metadata_conversion(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Update xarray dataset and data_vars with project-specific metadata fields.

    Parameters
    ----------
    d : xarray.Dataset
        Dataset with metadata to be updated.
    p : str
        Dataset project name.
    m : dict
        Metadata definition dictionary for project and variable(s).

    Returns
    -------
    xarray.Dataset
    """
    logging.info("Converting metadata to CF-like conventions.")

    header = m["Header"]

    # Static handling of version global attributes
    miranda_version = header.get("_miranda_version")
    if miranda_version:
        if isinstance(miranda_version, bool):
            header["miranda_version"] = __miranda_version__
        elif isinstance(miranda_version, dict):
            if p in miranda_version.keys():
                header["miranda_version"] = __miranda_version__
        else:
            logging.warning(
                f"`_miranda_version` not set for project `{p}`. Not appending."
            )
    if "_miranda_version" in header:
        del header["_miranda_version"]

    frequency = m["Header"].get("_frequency")
    if frequency:
        if isinstance(frequency, bool):
            _, m["Header"]["frequency"] = get_time_frequency(d)
        elif isinstance(frequency, dict):
            if p in frequency.keys():
                m["Header"]["frequency"] = get_time_frequency(d)
        else:
            logging.warning("`frequency` not set for project. Not appending.")
    if "_frequency" in m["Header"]:
        del m["Header"]["_frequency"]

    # Conditional handling of global attributes based on project name
    for field in [f for f in header if f.startswith("_")]:
        if isinstance(header[field], list):
            if p in header[field]:
                attr_treatment = header[field][p]
            else:
                logging.warning(
                    f"Attribute handling (`{field}`) not set for project `{p}`. Continuing..."
                )
                continue
        elif isinstance(header[field], dict):
            attr_treatment = header[field]
        else:
            raise AttributeError(
                "Attribute handling configuration seems to not be properly configured. Verify JSON."
            )

        if field == "_map_attrs":
            for attribute, mapping in attr_treatment.items():
                header[mapping] = d.attrs[attribute]
                del d.attrs[attribute]
        elif field == "_remove_attrs":
            for ff in attr_treatment:
                del d.attrs[ff]
        elif isinstance(attr_treatment, str):
            if field[1:] in d.attrs:
                logging.warning(
                    f"Overwriting `{field[1:]}` based on JSON configuration."
                )
            header[field[1:]] = attr_treatment
        else:
            raise AttributeError(
                f"Attribute treatment `{field}` is not properly configured. Verify JSON."
            )

        del header[field]

    # Add global attributes
    d.attrs.update(header)
    d.attrs.update(dict(project=p))

    # Date-based versioning
    if not d.attrs.get("version"):
        d.attrs.update(dict(version=f"v{VERSION}"))

    prev_history = d.attrs.get("history", "")
    history = (
        f"[{datetime.datetime.now()}] "
        "Converted variables and modified metadata for CF-like compliance: "
        f"{prev_history}".strip()
    )
    d.attrs.update(dict(history=history))

    return d


def dataset_corrections(ds: xr.Dataset, project: str) -> xr.Dataset:
    """Convert variables to CF-compliant format"""
    metadata_definition = load_json_data_mappings(project)

    ds = _correct_units_names(ds, project, metadata_definition)
    ds = _transform(ds, project, metadata_definition)
    ds = _invert_sign(ds, project, metadata_definition)
    ds = _units_cf_conversion(ds, metadata_definition)
    ds = _clip_values(ds, project, metadata_definition)

    ds = dims_conversion(ds, project, metadata_definition)
    ds = _ensure_correct_time(ds, project, metadata_definition)
    ds = _offset_time(ds, project, metadata_definition)

    ds = variable_conversion(ds, project, metadata_definition)

    ds = metadata_conversion(ds, project, metadata_definition)

    ds.attrs["history"] = (
        f"{datetime.datetime.now()}: "
        f"Variables converted from original files using miranda.convert.{dataset_corrections.__name__}. "
        f"{ds.attrs.get('history')}".strip()
    )

    return ds


def dataset_conversion(
    input_files: (
        str
        | os.PathLike
        | Sequence[str | os.PathLike]
        | Iterator[os.PathLike]
        | xr.Dataset
    ),
    project: str,
    domain: str | None = None,
    mask: xr.Dataset | xr.DataArray | None = None,
    mask_cutoff: float | bool = False,
    regrid: bool = False,
    add_version_hashes: bool = True,
    preprocess: Callable | str | None = "auto",
    **xr_kwargs,
) -> xr.Dataset | xr.DataArray:
    """Convert an existing Xarray-compatible dataset to another format with variable corrections applied.

    Parameters
    ----------
    input_files : str or os.PathLike or Sequence[str or os.PathLike] or Iterator[os.PathLike] or xr.Dataset
        Files or objects to be converted.
        If sent a list or GeneratorType, will open with :py:func:`xarray.open_mfdataset` and concatenate files.
    project : {"cordex", "cmip5", "cmip6", "ets-grnch", "isimip-ft", "pcic-candcs-u6", "converted"}
        Project name for decoding/handling purposes.
    domain: {"global", "nam", "can", "qc", "mtl"}, optional
        Domain to perform subsetting for. Default: None.
    mask : Optional[Union[xr.Dataset, xr.DataArray]]
        DataArray or single data_variable dataset containing mask.
    mask_cutoff : float or bool
        If land_sea_mask supplied, the threshold above which to mask with land_sea_mask. Default: False.
    regrid : bool
        Performing regridding with xesmf. Default: False.
    add_version_hashes : bool
        If True, version name and sha256sum of source file(s) will be added as a field among the global attributes.
    preprocess : callable or str, optional
        Preprocessing functions to perform over each Dataset.
        Default: "auto" - Run preprocessing fixes based on supplied fields from metadata definition.
        Callable - Runs function over Dataset (single) or supplied to `preprocess` (multifile dataset).
    **xr_kwargs
        Arguments passed directly to xarray.

    Returns
    -------
    xr.Dataset or xr.DataArray
    """
    if isinstance(input_files, xr.Dataset):
        ds = input_files
    else:
        if isinstance(input_files, (str, os.PathLike)):
            if Path(input_files).is_dir():
                files = []
                files.extend([f for f in Path(input_files).glob("*.nc")])
                files.extend([f for f in Path(input_files).glob("*.zarr")])
            else:
                files = [Path(input_files)]
        elif isinstance(input_files, (Sequence, Iterator)):
            files = [Path(f) for f in input_files]
        else:
            files = input_files
        version_hashes = dict()
        if add_version_hashes:
            for file in files:
                version_hashes[file.name] = find_version_hash(file)

        preprocess_kwargs = dict()
        if preprocess:
            if preprocess == "auto":
                preprocess_kwargs.update(
                    preprocess=partial(preprocessing_corrections, project=project)
                )
            elif isinstance(preprocess, Callable):
                preprocess_kwargs.update(preprocess=preprocess)

        if len(files) == 1:
            ds = xr.open_dataset(files[0], **xr_kwargs)
            for _, process in preprocess_kwargs.items():
                ds = process(ds)
        else:
            ds = xr.open_mfdataset(files, **xr_kwargs, **preprocess_kwargs)
        if version_hashes:
            ds.attrs.update(dict(original_files=str(version_hashes)))

    ds = dataset_corrections(ds, project)

    if domain:
        ds = subset_domain(ds, domain)

    if isinstance(mask, (str, Path)):
        mask = xr.open_dataset(mask)
    if isinstance(mask, (xr.Dataset, xr.DataArray)):
        if regrid:
            mask = conservative_regrid(ds, mask)
        ds = threshold_mask(ds, mask=mask, mask_cutoff=mask_cutoff)

    return ds
