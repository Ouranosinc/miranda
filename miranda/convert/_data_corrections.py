import datetime
import json
import logging.config
import os
import shutil
from pathlib import Path
from typing import Dict, Iterator, Optional, Sequence, Union

import numpy as np
import xarray as xr
from xclim.core import units

from miranda import __version__ as __miranda_version__
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

from .utils import delayed_write, find_version_hash

logging.config.dictConfig(LOGGING_CONFIG)

LATLON_COORDINATE_PRECISION = dict()
LATLON_COORDINATE_PRECISION["era5-land"] = 4

VERSION = datetime.datetime.now().strftime("%Y.%m.%d")

__all__ = ["file_conversion", "load_json_data_mappings", "variable_conversion"]


def load_json_data_mappings(project: str) -> dict:
    data_folder = Path(__file__).parent / "data"

    if project.startswith("era5"):
        metadata_definition = json.load(open(data_folder / "ecmwf_cf_attrs.json"))
    elif project in ["agcfsr", "agmerra2"]:  # This should handle the AG versions:
        metadata_definition = json.load(open(data_folder / "nasa_cf_attrs.json"))
    elif project in ["cordex", "cmip5", "cmip6"]:
        metadata_definition = json.load(open(data_folder / "cmip_ouranos_attrs.json"))
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
        if m["variable_entry"].get(v):
            if m["variable_entry"][v].get(key):
                if p in m["variable_entry"][v][key].keys():
                    print(m["variable_entry"][v][key][p])
                    d[v].attrs["units"] = m["variable_entry"][v][key][p]

    if m["variable_entry"].get("time"):
        if m["variable_entry"]["time"].get(key):
            if p in m["variable_entry"]["time"][key].keys():
                d["time"].attrs["units"] = m["variable_entry"]["time"][key][p]

    return d


# for de-accumulation or conversion to flux
def _transform(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    key = "_transformation"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    for vv in d.data_vars:
        converted = False
        if m["variable_entry"].get(vv):
            if m["variable_entry"][vv].get(key):
                if p in m["variable_entry"][vv][key].keys():
                    print("here!")

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
        if m["variable_entry"].get(vv):
            if m["variable_entry"][vv].get(key):
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
                            out["time"] = out.time - np.timedelta64(
                                offset[0], offset[1]
                            )
                            d_out[out.name] = out
                            converted = True
                    else:
                        logging.info(
                            f"No time offsetting needed for `{vv}` in `{p}` (Explicitly set to False)."
                        )
                else:
                    logging.info(
                        f"No time offsetting needed for `{vv}` in project `{p}`."
                    )
        if not converted:
            d_out[vv] = d[vv]
    return d_out


def _invert_sign(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    key = "_invert_sign"
    d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
    for vv in d.data_vars:
        converted = False
        if m["variable_entry"].get(vv):
            if m["variable_entry"][vv].get(key):
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
        if m["variable_entry"]["time"].get("units"):
            d["time"]["units"] = m["variable_entry"]["time"]["units"]

    for v in d.data_vars:
        if descriptions.get(v):
            d[v] = units.convert_units_to(d[v], descriptions[v]["units"])

    return d


# Add and update existing metadata fields
def metadata_conversion(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    logging.info("Converting metadata to CF-like conventions.")

    # Static handling of version global attributes
    miranda_version = m["Header"].get("_miranda_version")
    if miranda_version:
        if isinstance(miranda_version, bool):
            m["Header"]["miranda_version"] = __miranda_version__
        elif isinstance(miranda_version, dict):
            if p in miranda_version.keys():
                m["Header"]["miranda_version"] = __miranda_version__
        else:
            logging.warning("`__miranda_version__` not set for project. Not appending.")
    if "_miranda_version" in m["Header"]:
        del m["Header"]["_miranda_version"]

    # Conditional handling of global attributes based on project name
    cond_header = ["source", "doi"]
    for field in cond_header:
        if f"_{field}" in m["Header"].keys():
            if p in m["Header"][f"_{field}"].keys():
                m["Header"][field] = m["Header"][f"_{field}"][p]
            elif field in m["Header"].keys():
                pass
            else:
                raise AttributeError(f"`{field}` not found for project dataset.")
            del m["Header"][f"_{field}"]

    # Add global attributes
    d.attrs.update(m["Header"])
    d.attrs.update(dict(project=p))

    # Date-based versioning
    if not d.attrs.get("version"):
        d.attrs.update(dict(version=f"v{VERSION}"))

    prev_history = d.attrs.get("history", "")
    history = (
        f"[{datetime.datetime.now()}] "
        "Converted variables and modified metadata for CF-like compliance."
        f" {prev_history}".strip()
    )
    d.attrs.update(dict(history=history))
    descriptions = m["variable_entry"]

    time_correction_fields = ["_corrected_units", "_ensure_correct_time"]

    if "time" in m["variable_entry"].keys():
        for field in time_correction_fields:
            if field in m["variable_entry"]["time"].keys():
                del descriptions["time"][field]
        d["time"].attrs.update(descriptions["time"])

    # Add variable metadata and remove nonstandard entries
    data_vars_correction_fields = [
        "_corrected_units",
        "_invert_sign",
        "_offset_time",
        "_transformation",
    ]
    for v in d.data_vars:
        if v in descriptions.keys():
            for field in data_vars_correction_fields:
                if field in descriptions[v].keys():
                    del descriptions[v][field]
            d[v].attrs.update(descriptions[v])

    # Rename data variables
    for v in d.data_vars:
        if v in descriptions.keys():
            try:
                cf_name = descriptions[v]["_cf_variable_name"]
                d = d.rename({v: cf_name})
                d[cf_name].attrs.update(dict(original_variable=v))
                del d[cf_name].attrs["_cf_variable_name"]
            except (ValueError, IndexError):
                pass
    return d


def _ensure_correct_time(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    key = "_ensure_correct_time"

    if "time" not in m["variable_entry"].keys():
        logging.warning(f"No time corrections listed for project `{p}`. Continuing...")
        return d

    if "time" not in d.data_vars:
        logging.info(
            "No time dimension among data variables: "
            f"{' ,'.join([str(v) for v in d.data_vars])}. "
            "Continuing..."
        )
        return d

    if key in m["variable_entry"]["time"].keys():
        correct_times = m["variable_entry"]["time"][key][p]
        freq_found = xr.infer_freq(d.time)
        if not freq_found:
            raise ValueError(
                "Time frequency could not be found. There may be missing timesteps."
            )

        if freq_found in ["M", "A"]:
            freq_found = f"{freq_found}S"

        if freq_found not in correct_times:
            raise ValueError(
                f"Time frequency {freq_found} not among allowed frequencies: "
                f"{' ,'.join(correct_times)} for project `{p}`."
            )

        logging.info(f"Resampling dataset with time frequency: {freq_found}.")
        d_out = d.assign_coords(
            time=d.time.resample(time=freq_found).mean(dim="time").time
        )

        prev_history = d.attrs.get("history", "")
        history = f"[{datetime.datetime.now()}] Resampled time with `freq={freq_found}`.{prev_history}"
        d.attrs.update(dict(history=history))
        d = d_out

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
        except (KeyError, ValueError):
            pass
        if p in LATLON_COORDINATE_PRECISION.keys():
            d[new] = d[new].round(LATLON_COORDINATE_PRECISION[p])
    if sort_dims:
        d = d.sortby(sort_dims)
    return d


def variable_conversion(ds: xr.Dataset, project: str) -> xr.Dataset:
    """Convert variables to CF-compliant format"""
    metadata_definition = load_json_data_mappings(project)

    ds = _correct_units_names(ds, project, metadata_definition)
    ds = _transform(ds, project, metadata_definition)
    ds = _offset_time(ds, project, metadata_definition)
    ds = _invert_sign(ds, project, metadata_definition)
    ds = _units_cf_conversion(ds, metadata_definition)
    ds = _ensure_correct_time(ds, project, metadata_definition)
    ds = _dims_conversion(ds, project)

    ds = metadata_conversion(ds, project, metadata_definition)

    return ds


def file_conversion(
    files: Union[
        str, os.PathLike, Sequence[Union[str, os.PathLike]], Iterator[os.PathLike]
    ],
    project: str,
    output_path: Union[str, os.PathLike],
    output_format: str,
    chunks: Optional[dict] = None,
    overwrite: bool = False,
    add_version_hashes: bool = True,
    compute: bool = True,
    **xr_kwargs,
) -> None:
    """Convert an existing Xarray-compatible dataset to another format with variable corrections applied.

    Parameters
    ----------
    files : str or os.PathLike or Sequence[str or os.PathLike] or Iterator[os.PathLike]
        Files to be converted.
        If sent a list or GeneratorType, will open with :py:func:`xarray.open_mfdataset` and concatenate files.
    project : {"cordex", "cmip5", "cmip6", "isimip-ft", "pcic-candcs-u6", "converted"}

    output_path : str or os.PathLike
        Output folder path.
    output_format: {"netcdf", "zarr"}
        Output data container type.
    chunks : dict, optional
        Chunking layout to be written to new files. If None, chunking will be left to the relevant backend engine.
    overwrite: bool
        Whether to remove existing files or fail if files already exist.
    add_version_hashes: bool
        If True, version name and sha256sum of source file(s) will be added as a field among the global attributes.
    compute: bool
        If True, files will be converted with each call to file conversion.
        If False, will return a dask.Delayed object that can be computed later.
        Default: True.
    **xr_kwargs
        Arguments passed directly to xarray.

    Returns
    -------
    dask.Delayed or None
    """
    if output_format.lower() not in {"netcdf", "zarr"}:
        raise NotImplementedError(f"Format: {output_format}.")
    else:
        suffix = dict(netcdf="nc", zarr="zarr")[output_format]

    if isinstance(output_path, str):
        output_path = Path(output_path)

    if chunks is None:
        output_chunks = dict()
    else:
        output_chunks = chunks

    if isinstance(files, (str, os.PathLike)):
        files = [Path(files)]
    elif isinstance(files, (Sequence, Iterator)):
        files = [Path(f) for f in files]

    version_hashes = dict()
    if add_version_hashes:
        for file in files:
            version_hashes[file.name] = find_version_hash(file)

    if len(files) == 1:
        ds = xr.open_dataset(files[0], **xr_kwargs)
    else:
        ds = xr.open_mfdataset(files, **xr_kwargs)

    if version_hashes:
        ds.attrs.update(dict(original_files=str(version_hashes)))

    ds = variable_conversion(ds, project)
    outfile = f"{Path(files[0].stem)}.{suffix}"
    outfile_path = output_path.joinpath(outfile)

    if overwrite and outfile_path.exists():
        logging.warning(f"Removing existing {output_format} files for {files[0].name}.")
        if outfile_path.is_dir():
            shutil.rmtree(outfile_path)
        if outfile_path.is_file():
            outfile_path.unlink()

    write_object = delayed_write(
        ds, outfile_path, output_chunks, output_format, overwrite
    )
    if compute:
        return write_object.compute()
    return write_object
