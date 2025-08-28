from __future__ import annotations
import logging

import xarray as xr
import xclim.core.units
from xclim.core import units

from miranda.treatments.utils import (
    _get_section_entry_key,  # noqa
    _iter_entry_key,  # noqa
)
from miranda.units import check_time_frequency


logger = logging.getLogger("miranda.treatments.variables")

__all__ = [
    "cf_units_conversion",
    "clip_values",
    "correct_unit_names",
    "invert_value_sign",
    "transform_values",
    "variable_conversion",
]


def correct_unit_names(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Correct unit names."""
    key = "_corrected_units"
    for var, val in _iter_entry_key(d, m, "variables", key, p):
        if val:
            d[var].attrs["units"] = val
            prev_history = d.attrs.get("history", "")
            history = f"Corrected units name for variable `{var}` to `{val}`. {prev_history}"
            d.attrs.update(dict(history=history))

    return d


# for de-accumulation or conversion to flux
def transform_values(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Transform dataset values according to operation listed."""
    key = "_transformation"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    converted = []
    offset, offset_meaning = None, None

    time_freq = dict()
    expected_period = _get_section_entry_key(m, "dimensions", "time", "_ensure_correct_time", p)
    if isinstance(expected_period, str):
        time_freq["expected_period"] = expected_period

    for vv, trans in _iter_entry_key(d, m, "variables", key, p):
        if trans:
            if trans == "deaccumulate":
                # Time-step accumulated total to time-based flux (de-accumulation)
                if offset is None and offset_meaning is None:
                    try:
                        offset, offset_meaning = check_time_frequency(d, **time_freq)
                    except TypeError:
                        logger.error("Unable to parse the time frequency. Verify data integrity before retrying.")
                        raise

                msg = f"De-accumulating units for variable `{vv}`."
                logger.info(msg)
                with xr.set_options(keep_attrs=True):
                    out = d[vv].diff(dim="time")
                    out = d[vv].where(
                        getattr(d[vv].time.dt, offset_meaning) == offset[0],
                        out.broadcast_like(d[vv]),
                    )
                    out = units.amount2rate(out)
                    d_out[vv] = out
                converted.append(vv)
            elif trans == "amount2rate":
                # NOTE: This treatment is no longer needed in xclim v0.43.0+ but is kept for backwards compatibility
                # frequency-based totals to time-based flux
                msg = f"Performing amount-to-rate units conversion for variable `{vv}`."
                logger.info(msg)
                with xr.set_options(keep_attrs=True):
                    out = units.amount2rate(d[vv])
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
                            raise NotImplementedError(f"Op transform doesn't implement the «{op}» operator.")
                converted.append(vv)
            else:
                raise NotImplementedError(f"Unknown transformation: {trans}")
        elif trans is False:
            msg = f"No transformations needed for `{vv}` (Explicitly set to False)."
            logger.info(msg)
            continue

        prev_history = d.attrs.get("history", "")
        history = f"Transformed variable `{vv}` values using method `{trans}`. {prev_history}"
        d_out.attrs.update(dict(history=history))

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]
    return d_out


def invert_value_sign(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Flip value of DataArray."""
    key = "_invert_sign"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    converted = []
    for vv, inv_sign in _iter_entry_key(d, m, "variables", key, p):
        if inv_sign:
            msg = f"Inverting sign for `{vv}` (switching direction of values)."
            logger.info(msg)
            with xr.set_options(keep_attrs=True):
                out = d[vv]
                d_out[out.name] = -out
            converted.append(vv)
        elif inv_sign is False:
            msg = f"No sign inversion needed for `{vv}` in `{p}` (Explicitly set to False)."
            logger.info(msg)
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
def cf_units_conversion(d: xr.Dataset, m: dict) -> xr.Dataset:
    """Perform pint-based units-conversion."""
    if "time" in m["dimensions"].keys():
        if m["dimensions"]["time"].get("units"):
            d["time"]["units"] = m["dimensions"]["time"]["units"]

    for vv, unit in _iter_entry_key(d, m, "variables", "units", None):
        context = m["variables"][vv].get("_units_context", None)
        if unit:
            with xr.set_options(keep_attrs=True):
                d[vv] = units.convert_units_to(d[vv], unit, context=context)
            prev_history = d.attrs.get("history", "")
            history = f"Converted variable `{vv}` to CF-compliant units (`{unit}`). {prev_history}"
            d.attrs.update(dict(history=history))

    return d


# For clipping variable values to an established maximum/minimum
def clip_values(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Clip values to an appropriate range,."""
    key = "_clip_values"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    converted = []
    for vv in d.data_vars:
        if vv in m["variables"].keys():
            clip_vals = _get_section_entry_key(m, "variables", vv, key, p)
            if clip_vals:
                min_value, max_value = None, None
                # Gather unit conversion context, if applicable
                context = clip_vals.get("context", None)
                for op, value in clip_vals.items():
                    if op == "min":
                        min_value = xclim.core.units.convert_units_to(value, d[vv], context)
                    if op == "max":
                        max_value = xclim.core.units.convert_units_to(value, d[vv], context)
                msg = f"Clipping min/max values for `{vv}` ({min_value}/{max_value})."
                logger.info(msg)
                with xr.set_options(keep_attrs=True):
                    out = d[vv]
                    d_out[out.name] = out.clip(min_value, max_value)
                converted.append(vv)
            elif clip_values is False:
                msg = f"No clipping of values needed for `{vv}` in `{p}` (Explicitly set to False)."
                logger.info(msg)
                continue
            else:
                msg = f"Unknown clipping values for `{vv}` in `{p}`."
                logger.info(msg)
                continue

            prev_history = d.attrs.get("history", "")
            history = f"Clipped variable `{vv}` with `min={min_value}` and `max={max_value}`. {prev_history}"
            d_out.attrs.update(dict(history=history))

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]

    return d_out


# For renaming and reordering lat and lon dims


def variable_conversion(d: xr.Dataset, p: str | None, m: dict) -> xr.Dataset:
    """
    Add variable metadata and remove nonstandard entries.

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
    for orig_var_name, cf_name in _iter_entry_key(d, m, "variables", "_cf_variable_name", p):
        if cf_name is not None:
            d = d.rename({orig_var_name: cf_name})
            d[cf_name].attrs.update(dict(original_variable=orig_var_name))
            del d[cf_name].attrs["_cf_variable_name"]

    return d
