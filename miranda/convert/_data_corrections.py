import datetime
import json
import logging.config
import os
import shutil
from functools import partial
from pathlib import Path
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Union

import numpy as np
import xarray as xr
from clisops.core.subset import subset_bbox
from xarray.coding import times
from xclim.core import units
from xclim.core.calendar import parse_offset

from miranda import __version__ as __miranda_version__
from miranda.decode import date_parser
from miranda.gis import subset_domain
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

from .utils import delayed_write, find_version_hash, name_output_file

logging.config.dictConfig(LOGGING_CONFIG)

VERSION = datetime.datetime.now().strftime("%Y.%m.%d")

__all__ = [
    "dataset_corrections",
    "dims_conversion",
    "file_conversion",
    "load_json_data_mappings",
    "metadata_conversion",
    "threshold_land_sea_mask",
    "variable_conversion",
]


def load_json_data_mappings(project: str) -> dict:
    data_folder = Path(__file__).parent / "data"

    if project.startswith("era5"):
        metadata_definition = json.load(open(data_folder / "ecmwf_cf_attrs.json"))
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
        metadata_definition = json.load(open(data_folder / "ec_cf_attrs.json"))
    else:
        raise NotImplementedError()

    return metadata_definition


def _get_section_entry_key(meta, entry, var, key, project):
    var_meta = meta[entry].get(var, {})
    if key in var_meta:
        if isinstance(var_meta[key], dict):
            return var_meta[key].get(project)
        return var_meta[key]
    return None


def _iter_entry_key(ds, meta, entry, key, project):
    for vv in set(ds.data_vars).intersection(meta[entry]):
        val = _get_section_entry_key(meta, entry, vv, key, project)
        yield vv, val


def _simple_fix_dims(
    d: Union[xr.Dataset, xr.DataArray]
) -> Union[xr.Dataset, xr.DataArray]:
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
    ds: Union[xr.DataArray, xr.Dataset], ref_grid: Union[xr.DataArray, xr.Dataset]
) -> Union[xr.DataArray, xr.Dataset]:
    """Perform a conservative_normed regridding"""
    try:
        import xesmf as xe
    except ModuleNotFoundError:
        logging.warning(
            "This function requires the `xesmf` library which is not installed. "
            "Regridding step will be skipped."
        )
        return ds

    ref_grid = _simple_fix_dims(ref_grid)
    method = "conservative_normed"

    regridder = xe.Regridder(ds, ref_grid, method, periodic=False)
    ds = regridder(ds)

    ds.attrs["history"] = (
        f"{datetime.datetime.now()}:"
        f"Regridded dataset using xesmf with method: {method}. "
        f"{ds.attrs.get('history')}".strip()
    )
    return ds


def threshold_land_sea_mask(
    ds: Union[xr.Dataset, xr.DataArray],
    *,
    land_sea_mask: Union[xr.Dataset, xr.DataArray],
    land_sea_cutoff: float = 0.5,
) -> Union[xr.Dataset, xr.DataArray]:
    """Land-Sea mask operations.

    Parameters
    ----------
    ds : Union[xr.Dataset, str, os.PathLike]
    land_sea_mask : Union[xr.Dataset, xr.DataArray]
    land_sea_cutoff : float

    Returns
    -------
    Union[xr.Dataset, xr.DataArray]
    """
    logging.info(
        f"Masking variable with land-sea mask at `{land_sea_cutoff}` cutoff value."
    )

    land_sea_mask = _simple_fix_dims(land_sea_mask)

    if isinstance(land_sea_mask, xr.Dataset):
        if len(land_sea_mask.data_vars) == 1:
            land_sea_mask = land_sea_mask[list(land_sea_mask.data_vars)[0]]
        else:
            raise ValueError(
                "More than one data variable found in land-sea mask. Supply a DataArray instead."
            )

    lon_bounds = np.array([ds.lon.min(), ds.lon.max()])
    lat_bounds = np.array([ds.lat.min(), ds.lat.max()])

    lsm = subset_bbox(
        land_sea_mask,
        lon_bnds=lon_bounds,
        lat_bnds=lat_bounds,
    ).load()

    lsm = lsm.where(land_sea_mask >= land_sea_cutoff)
    ds = ds.where(lsm.notnull())

    if lsm.min() >= 0:
        if lsm.max() <= 1.00000001:
            ds.attrs["land_sea_cutoff"] = f"{land_sea_cutoff * 100} %"
        elif lsm.max() <= 100.00000001:
            ds.attrs["land_sea_cutoff"] = f"{land_sea_cutoff} %"
    else:
        ds.attrs["land_sea_cutoff"] = f"{land_sea_cutoff}"

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
    time = xr.coding.times.decode_cf_datetime(
        vals, units=f"days since {date}", calendar="standard"
    )
    d = d.assign_coords({field: time})
    return d


def correct_var_names(d: xr.Dataset, split: str = "_", location: int = 0) -> xr.Dataset:
    filename = d.encoding["source"]
    new_name = Path(filename).stem.split(split)[location]
    old_name = list(d.data_vars.keys())[0]
    return d.rename({old_name: new_name})


def preprocess_corrections(ds: xr.Dataset, *, project: str) -> xr.Dataset:
    def _preprocess_correct(d: xr.Dataset, *, ops: List[partial]) -> xr.Dataset:
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


def _correct_units_names(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
    key = "_corrected_units"
    for var, val in _iter_entry_key(d, m, "variables", key, p):
        if val:
            d[var].attrs["units"] = val

    val_time = _get_section_entry_key(m, "variables", "time", key, p)
    if val_time:
        d["time"].attrs["units"] = val_time

    return d


# for de-accumulation or conversion to flux
def _transform(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
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
    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]
    return d_out


def _offset_time(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
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

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]
    return d_out


def _invert_sign(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
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

    # Copy unconverted variables
    for vv in d.data_vars:
        if vv not in converted:
            d_out[vv] = d[vv]
    return d_out


# For converting variable units to standard workflow units
def _units_cf_conversion(d: xr.Dataset, m: Dict) -> xr.Dataset:
    if "time" in m["dimensions"].keys():
        if m["dimensions"]["time"].get("units"):
            d["time"]["units"] = m["dimensions"]["time"]["units"]

    for vv, uni in _iter_entry_key(d, m, "variables", "units", None):
        if uni:
            with xr.set_options(keep_attrs=True):
                d[vv] = units.convert_units_to(d[vv], uni)

    return d


def _ensure_correct_time(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
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

        if freq_found in ["M", "A"]:
            freq_found = f"{freq_found}S"

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
                f"{' ,'.join(correct_times)}"
            )
            if isinstance(correct_time_entry, dict):
                error_msg = f"{error_msg} for project `{p}`."
            else:
                error_msg = f"{error_msg}."
            raise ValueError(error_msg)

        logging.info(f"Resampling dataset with time frequency: {freq_found}.")
        d_out = d.assign_coords(
            time=d.time.resample(time=freq_found).mean(dim="time").time
        )

        prev_history = d.attrs.get("history", "")
        history = f"[{datetime.datetime.now()}] Resampled time with `freq={freq_found}`.{prev_history}"
        d.attrs.update(dict(history=history))
        d = d_out

    return d


# For renaming and reordering lat and lon dims
def dims_conversion(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    # Rename dimensions to CF to their equivalents
    rename_dims = dict()
    for dim in d.dims:
        if dim in m["dimensions"].keys():
            cf_name = _get_section_entry_key(
                m, "dimensions", dim, "_cf_dimension_name", p
            )
            if cf_name:
                rename_dims[dim] = cf_name
    d = d.rename(rename_dims)

    # Perform lon wrapping if needed and sort dimensions
    sort_dims = []
    for new in ["lon", "lat"]:
        if new == "lon" and "lon" in d.dims:
            if np.any(d.lon > 180):
                lon1 = d.lon.where(d.lon <= 180.0, d.lon - 360.0)
                d[new] = lon1
        sort_dims.append(new)
        coord_precision = _get_section_entry_key(m, "dimensions", new, "_precision", p)
        if coord_precision is not None:
            d[new] = d[new].round(coord_precision)
    if sort_dims:
        d = d.sortby(sort_dims)

    # Ensure that lon and lat are written in proper order for plotting purposes
    transpose_order = ["lat", "lon"]
    if "time" in d.dims:
        transpose_order.insert(0, "time")
    transpose_order.extend(list(set(d.dims) - set(transpose_order)))
    d = d.transpose(*transpose_order)

    # Add dimension original name and update attrs
    dim_descriptions = m["dimensions"]
    for dim in m["dimensions"].keys():
        cf_name = dim_descriptions[dim].get("_cf_dimension_name")
        if cf_name is not None and cf_name in d.dims:
            d[cf_name].attrs.update(dict(original_variable=dim))
            for field in dim_descriptions[dim].keys():
                if not field.startswith("_"):
                    d[cf_name].attrs.update({field: dim_descriptions[dim][field]})

    return d


def variable_conversion(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    # Add variable metadata and remove nonstandard entries
    var_descriptions = m["variables"]
    var_correction_fields = [
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


def metadata_conversion(d: xr.Dataset, p: str, m: Dict) -> xr.Dataset:
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
    for field in [f for f in m["Header"] if f.startswith("_")]:
        if p in m["Header"][field]:
            m["Header"][field[1:]] = m["Header"][field][p]
        elif field[1:] in m["Header"]:
            pass
        else:
            raise AttributeError(f"`{field[1:]}` not found for project dataset.")
        del m["Header"][field]

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

    return d


def dataset_corrections(ds: xr.Dataset, project: str) -> xr.Dataset:
    """Convert variables to CF-compliant format"""
    metadata_definition = load_json_data_mappings(project)

    ds = _correct_units_names(ds, project, metadata_definition)
    ds = _transform(ds, project, metadata_definition)
    ds = _invert_sign(ds, project, metadata_definition)
    ds = _units_cf_conversion(ds, metadata_definition)

    ds = dims_conversion(ds, project, metadata_definition)
    ds = _ensure_correct_time(ds, project, metadata_definition)
    ds = _offset_time(ds, project, metadata_definition)

    ds = variable_conversion(ds, project, metadata_definition)

    ds = metadata_conversion(ds, project, metadata_definition)

    return ds


def file_conversion(
    files: Union[
        str, os.PathLike, Sequence[Union[str, os.PathLike]], Iterator[os.PathLike]
    ],
    project: str,
    output_path: Union[str, os.PathLike],
    output_format: str,
    *,
    domain: Optional[str] = None,
    land_sea_mask: Optional[Union[xr.Dataset, xr.DataArray]] = None,
    land_sea_cutoff: float = 0.5,
    chunks: Optional[dict] = None,
    overwrite: bool = False,
    add_version_hashes: bool = True,
    preprocess: Optional[Union[Callable, str]] = "auto",
    compute: bool = True,
    **xr_kwargs,
) -> Dict:
    """Convert an existing Xarray-compatible dataset to another format with variable corrections applied.

    Parameters
    ----------
    files : str or os.PathLike or Sequence[str or os.PathLike] or Iterator[os.PathLike]
        Files to be converted.
        If sent a list or GeneratorType, will open with :py:func:`xarray.open_mfdataset` and concatenate files.
    project : {"cordex", "cmip5", "cmip6", "ets-grnch", "isimip-ft", "pcic-candcs-u6", "converted"}
        Project name for decoding/handling purposes.
    output_path : str or os.PathLike
        Output folder path.
    output_format: {"netcdf", "zarr"}
        Output data container type.
    domain: {"global", "nam", "can", "qc", "mtl"}, optional
        Domain to perform subsetting for. Default: None.
    land_sea_mask : Optional[Union[xr.Dataset, xr.DataArray]]
        DataArray or single data_variable dataset containing land-sea mask.
    land_sea_cutoff : float
        If land_sea_mask supplied, the threshold above which to mask with land_sea_mask. Default: 0.5.
    chunks : dict, optional
        Chunking layout to be written to new files. If None, chunking will be left to the relevant backend engine.
    overwrite : bool
        Whether to remove existing files or fail if files already exist.
    add_version_hashes : bool
        If True, version name and sha256sum of source file(s) will be added as a field among the global attributes.
    preprocess : callable or str, optional
        Preprocessing functions to perform over each Dataset.
        Default: "auto" - Run preprocessing fixes based on supplied fields from metadata definition.
        Callable - Runs function over Dataset (single) or supplied to `preprocess` (multifile dataset).
    compute : bool
        If True, files will be converted with each call to file conversion.
        If False, will return a dask.Delayed object that can be computed later.
        Default: True.
    **xr_kwargs
        Arguments passed directly to xarray.

    Returns
    -------
    dict
    """
    if isinstance(output_path, str):
        output_path = Path(output_path)

    if isinstance(files, (str, os.PathLike)):
        files = [Path(files)]
    elif isinstance(files, (Sequence, Iterator)):
        files = [Path(f) for f in files]

    version_hashes = dict()
    if add_version_hashes:
        for file in files:
            version_hashes[file.name] = find_version_hash(file)

    preprocess_kwargs = dict()
    if preprocess:
        if preprocess == "auto":
            preprocess_kwargs.update(
                preprocess=partial(preprocess_corrections, project=project)
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
    ds.attrs["history"] = (
        f"{datetime.datetime.now()}: "
        f"Variables converted from original files using miranda.convert.{file_conversion.__name__}. "
        f"{ds.attrs.get('history')}".strip()
    )

    if domain:
        ds = subset_domain(ds, domain)

    if isinstance(land_sea_mask, (xr.Dataset, xr.DataArray)):
        logging.info(
            "Land-sea mask supplied. Performing conservative-normed regridding and masking."
        )
        land_sea_mask = conservative_regrid(land_sea_mask, ds)
        ds = threshold_land_sea_mask(
            ds, land_sea_mask=land_sea_mask, land_sea_cutoff=land_sea_cutoff
        )

    outfile = name_output_file(ds, project, output_format)
    outfile_path = output_path.joinpath(outfile)

    if overwrite and outfile_path.exists():
        logging.warning(f"Removing existing {output_format} files for {outfile}.")
        if outfile_path.is_dir():
            shutil.rmtree(outfile_path)
        if outfile_path.is_file():
            outfile_path.unlink()

    logging.info(f"Writing {outfile}.")
    write_object = delayed_write(
        ds,
        outfile_path,
        output_format,
        overwrite,
        target_chunks=chunks,
    )
    if compute:
        write_object.compute()
        return dict(path=outfile_path)
    return dict(path=outfile_path, object=write_object)
