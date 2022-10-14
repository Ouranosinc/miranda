import datetime
import json
import logging.config
from pathlib import Path
from typing import Dict, Union

import numpy as np
import xarray as xr
from xclim.core import units

from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

logging.config.dictConfig(LOGGING_CONFIG)

LATLON_COORDINATE_PRECISION = dict()
LATLON_COORDINATE_PRECISION["era5-land"] = 4

VERSION = datetime.datetime.now().strftime("%Y.%m.%d")

__all__ = ["load_json_data_mappings", "variable_conversion"]


def load_json_data_mappings(project: str) -> dict:
    data_folder = Path(__file__).parent / "data"

    if project.startswith("era5"):
        metadata_definition = json.load(open(data_folder / "ecmwf_cf_attrs.json"))
    elif project in ["agcfsr", "agmerra2"]:  # This should handle the AG versions:
        metadata_definition = json.load(open(data_folder / "nasa_cf_attrs.json"))
    elif project == "nrcan-gridded-10km":
        raise NotImplementedError()
    elif project == "wfdei-gem-capa":
        metadata_definition = json.load(open(data_folder / "usask_cf_attrs.json"))
    else:
        raise NotImplementedError()

    return metadata_definition


def _correct_units_names(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    key = "_corrected_units"
    for v in d.data_vars:
        if p in m["variable_entry"][v][key].keys():
            d[v].attrs["units"] = m["variable_entry"][v][key][p]

    if "time" in m["variable_entry"].keys():
        if p in m["variable_entry"]["time"][key].keys():
            d["time"].attrs["units"] = m["variable_entry"]["time"][key][p]

    return d


# for de-accumulation or conversion to flux
def _transform(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    key = "_transformation"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    for vv in d.data_vars:
        converted = False
        if hasattr(m["variable_entry"][vv], key):
            if p in m["variable_entry"][vv][key].keys():
                try:
                    offset, offset_meaning = get_time_frequency(d)
                except TypeError:
                    logging.error(
                        f"Unable to parse the time frequency for variable `{vv}`. "
                        "Verify data integrity before retrying."
                    )
                    raise

                if m["variable_entry"][vv][key][p] == "deaccumulate":
                    # Time-step accumulated total to time-based flux (de-accumulation)
                    logging.info(f"De-accumulating units for variable `{vv}`.")
                    with xr.set_options(keep_attrs=True):
                        out = d[vv].diff(dim="time")
                        out = d[vv].where(
                            getattr(d[vv].time.dt, offset_meaning) == offset[0],
                            out.broadcast_like(d[vv]),
                        )
                        out = units.amount2rate(out)
                    d_out[out.name] = out
                    converted = True
                elif m["variable_entry"][vv][key][p] == "amount2rate":
                    # frequency-based totals to time-based flux
                    logging.info(
                        f"Performing amount-to-rate units conversion for variable `{vv}`."
                    )
                    out = units.amount2rate(
                        d[vv],
                        out_units=m["variable_entry"][vv]["units"],
                    )
                    d_out[out.name] = out
                    converted = True
                else:
                    raise NotImplementedError(
                        f"Unknown transformation: {m['variable_entry'][vv][key][p]}"
                    )
        if not converted:
            d_out[vv] = d[vv]
    return d_out


def _offset_time(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    key = "_offset_time"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    for vv in d.data_vars:
        converted = False
        if hasattr(m["variable_entry"][vv], key):
            if p in m["variable_entry"][vv][key].keys():
                try:
                    offset, offset_meaning = get_time_frequency(d)
                except TypeError:
                    logging.error(
                        f"Unable to parse the time frequency for variable `{vv}`. "
                        "Verify data integrity before retrying."
                    )
                    raise

                if m["variable_entry"][vv][key][p]:
                    # Offset time by value of one time-step
                    logging.info(
                        f"Offsetting data for `{vv}` by `{offset[0]} {offset_meaning}(s)`."
                    )
                    with xr.set_options(keep_attrs=True):
                        out = d[vv]
                        out["time"] = out.time - np.timedelta64(offset[0], offset[1])
                        d_out[out.name] = out
                        converted = True
                else:
                    logging.info(
                        f"No time offsetting needed for `{vv}` in `{p}` (Explicitly set to False)."
                    )
            else:
                logging.info(f"No time offsetting needed for `{vv}` in project `{p}`.")
        if not converted:
            d_out[vv] = d[vv]
    return d_out


def _invert_sign(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    key = "_invert_sign"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    for vv in d.data_vars:
        converted = False
        if hasattr(m["variable_entry"][vv], key):
            if p in m["variable_entry"][vv][key].keys():
                if m["variable_entry"][vv][key][p]:
                    logging.info(
                        f"Inverting sign for `{vv}` (switching direction of values)."
                    )
                    with xr.set_options(keep_attrs=True):
                        out = d[vv]
                        d_out[out.name] = out.__invert__()
                        converted = True
                else:
                    logging.info(
                        f"No sign inversion needed for `{vv}` in `{p}` (Explicitly set to False)."
                    )
            else:
                logging.info(f"No sign inversion needed for `{vv}` in `{p}`.")
        if not converted:
            d_out[vv] = d[vv]
    return d_out


# For converting variable units to standard workflow units
def _units_cf_conversion(d: xr.Dataset, m: Dict) -> xr.Dataset:
    descriptions = m["variable_entry"]

    if "time" in m["variable_entry"].keys():
        d["time"]["units"] = m["variable_entry"]["time"]["units"]

    for v in d.data_vars:
        d[v] = units.convert_units_to(d[v], descriptions[v]["units"])

    return d


# Add and update existing metadata fields
def _metadata_conversion(d: xr.Dataset, p: str, o: str, m: Dict) -> xr.Dataset:
    logging.info("Converting metadata to CF-like conventions.")

    # Conditional handling of source based on project name
    if "_source" in m["Header"].keys():
        if p in m["Header"]["_source"].keys():
            m["Header"]["source"] = m["Header"]["_source"][p]
        elif "source" in m["Header"].keys():
            pass
        else:
            raise AttributeError("Source not found for project dataset.")
        del m["Header"]["_source"]

    # Conditional handling of DOI based on project name
    if "_doi" in m["Header"].keys():
        if p in m["Header"]["_doi"].keys():
            m["Header"]["doi"] = m["Header"]["_doi"][p]
        elif "doi" in m["Header"].keys():
            pass
        else:
            logging.warning("DOI not found for project dataset. Skipping.")
        del m["Header"]["_doi"]

    # Add global attributes
    d.attrs.update(m["Header"])
    d.attrs.update(dict(project=p, format=o))

    # Date-based versioning
    d.attrs.update(dict(version=f"v{VERSION}"))

    if hasattr(d.attrs, "history"):
        prev_history = f" {d.attrs['history']}"
    else:
        prev_history = ""

    history = (
        f"[{datetime.datetime.now()}] Converted from original data to {o} "
        f"with modified metadata for CF-like compliance.{prev_history}"
    )
    d.attrs.update(dict(history=history))
    descriptions = m["variable_entry"]

    if "time" in m["variable_entry"].keys():
        del descriptions["time"]["_corrected_units"]

    # Add variable metadata and remove nonstandard entries
    correction_fields = [
        "_corrected_units",
        "_invert_sign",
        "_offset_time",
        "_transformation",
    ]
    for v in d.data_vars:
        for field in correction_fields:
            if field in descriptions[v].keys():
                del descriptions[v][field]
        d[v].attrs.update(descriptions[v])

    # Rename data variables
    for v in d.data_vars:
        try:
            cf_name = descriptions[v]["_cf_variable_name"]
            d = d.rename({v: cf_name})
            d[cf_name].attrs.update(dict(original_variable=v))
            del d[cf_name].attrs["_cf_variable_name"]
        except (ValueError, IndexError):
            pass
    return d


# For renaming lat and lon dims
def _dims_conversion(d: xr.Dataset, p: str) -> xr.Dataset:
    sort_dims = []
    for orig, new in dict(longitude="lon", latitude="lat").items():
        try:

            d = d.rename({orig: new})
            if new == "lon" and np.any(d.lon > 180):
                lon1 = d.lon.where(d.lon <= 180.0, d.lon - 360.0)
                d[new] = lon1
            sort_dims.append(new)
        except KeyError:
            pass
        if p in LATLON_COORDINATE_PRECISION.keys():
            d[new] = d[new].round(LATLON_COORDINATE_PRECISION[p])
    if sort_dims:
        d = d.sortby(sort_dims)
    return d


def variable_conversion(ds: xr.Dataset, project: str, output_format: str) -> xr.Dataset:
    """Convert variables to CF-compliant format"""
    metadata_definition = load_json_data_mappings(project)
    ds = _correct_units_names(ds, project, metadata_definition)
    ds = _transform(ds, project, metadata_definition)
    ds = _offset_time(ds, project, metadata_definition)
    ds = _invert_sign(ds, project, metadata_definition)
    ds = _units_cf_conversion(ds, metadata_definition)
    ds = _metadata_conversion(ds, project, output_format, metadata_definition)
    ds = _dims_conversion(ds, project)

    return ds
