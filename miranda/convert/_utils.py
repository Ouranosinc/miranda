import datetime
import json
import logging.config
import os
from pathlib import Path
from typing import Dict, Optional, Union

import netCDF4
import numpy as np
import xarray as xr
import zarr
from clisops.core import subset
from xclim.core import units
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "get_chunks_on_disk",
    "add_ar6_regions",
    "daily_aggregation",
    "delayed_write",
    "variable_conversion",
]

LATLON_COORDINATE_PRECISION = dict()
LATLON_COORDINATE_PRECISION["era5-land"] = 4

VERSION = datetime.datetime.now().strftime("%Y.%m.%d")


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


def get_chunks_on_disk(nc_file: Union[os.PathLike, str]) -> dict:
    """

    Parameters
    ----------
    nc_file: Path or str

    Returns
    -------
    dict
    """
    # FIXME: This does not support zarr
    # TODO: This needs to be optimized for dask. See: https://github.com/Ouranosinc/miranda/pull/24/files#r840617216
    ds = netCDF4.Dataset(nc_file)
    chunks = dict()
    for v in ds.variables:
        chunks[v] = dict()
        for ii, dim in enumerate(ds[v].dimensions):
            chunks[v][dim] = ds[v].chunking()[ii]
    return chunks


def add_ar6_regions(ds: xr.Dataset) -> xr.Dataset:
    """Add the IPCC AR6 Regions to dataset.

    Parameters
    ----------
    ds : xarray.Dataset

    Returns
    -------
    xarray.Dataset
    """
    try:
        import regionmask  # noqa
    except ImportError:
        raise ImportError(
            f"{add_ar6_regions.__name__} functions require additional dependencies. "
            "Please install them with `pip install miranda[full]`."
        )

    mask = regionmask.defined_regions.ar6.all.mask(ds.lon, ds.lat)
    ds = ds.assign_coords(region=mask)
    return ds


def variable_conversion(ds: xr.Dataset, project: str, output_format: str) -> xr.Dataset:
    """Convert variables to CF-compliant format"""

    def _correct_units_names(d: xr.Dataset, p: str, m: Dict):
        key = "_corrected_units"
        for v in d.data_vars:
            if p in m["variable_entry"][v][key].keys():
                d[v].attrs["units"] = m["variable_entry"][v][key][project]

        if "time" in m["variable_entry"].keys():
            if p in m["variable_entry"]["time"][key].keys():
                d["time"].attrs["units"] = m["variable_entry"]["time"][key][project]

        return d

    # for de-accumulation or conversion to flux
    def _transform(d: xr.Dataset, p: str, m: Dict):
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
                raise AttributeError("Source not found for project dataset")
            del m["Header"]["_source"]

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
                if hasattr(descriptions[v], field):
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
    def _dims_conversion(d: xr.Dataset):
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
            if project in LATLON_COORDINATE_PRECISION.keys():
                d[new] = d[new].round(LATLON_COORDINATE_PRECISION[project])
        if sort_dims:
            d = d.sortby(sort_dims)
        return d

    metadata_definition = load_json_data_mappings(project)
    ds = _correct_units_names(ds, project, metadata_definition)
    ds = _transform(ds, project, metadata_definition)
    ds = _offset_time(ds, project, metadata_definition)
    ds = _invert_sign(ds, project, metadata_definition)
    ds = _units_cf_conversion(ds, metadata_definition)
    ds = _metadata_conversion(ds, project, output_format, metadata_definition)
    ds = _dims_conversion(ds)

    return ds


def daily_aggregation(ds) -> Dict[str, xr.Dataset]:
    logging.info("Creating daily upscaled climate variables.")

    daily_dataset = dict()
    for variable in ds.data_vars:
        if variable in ["tas", "tdps"]:
            # Some looping to deal with memory consumption issues
            for v, func in {
                f"{variable}max": "max",
                f"{variable}min": "min",
                f"{variable}": "mean",
            }.items():
                ds_out = xr.Dataset()
                ds_out.attrs = ds.attrs.copy()
                ds_out.attrs["frequency"] = "day"

                method = (
                    f"time: {func}{'imum' if func != 'mean' else ''} (interval: 1 day)"
                )
                ds_out.attrs["cell_methods"] = method

                if v == "tas" and not hasattr(ds, "tas"):
                    ds_out[v] = tas(tasmax=ds.tasmax, tasmin=ds.tasmin)
                else:
                    # Thanks for the help, xclim contributors
                    r = ds[variable].resample(time="D")
                    ds_out[v] = getattr(r, func)(dim="time", keep_attrs=True)

                daily_dataset[v] = ds_out
                del ds_out

        elif variable in [
            "evspsblpot",
            "hfls",
            "hfss",
            "pr",
            "prsn",
            "rsds",
            "rlds",
            "snd",
            "snr",
            "snw",
        ]:
            ds_out = xr.Dataset()
            ds_out.attrs = ds.attrs.copy()
            ds_out.attrs["frequency"] = "day"
            ds_out.attrs["cell_methods"] = "time: mean (interval: 1 day)"
            logging.info(f"Converting {variable} to daily time step (daily mean).")
            ds_out[variable] = (
                ds[variable].resample(time="D").mean(dim="time", keep_attrs=True)
            )

            daily_dataset[variable] = ds_out
            del ds_out
        else:
            continue

    return daily_dataset


def threshold_land_sea_mask(
    ds: Union[xr.Dataset, str, os.PathLike],
    *,
    land_sea_mask: Dict[str, Union[os.PathLike, str]],
    land_sea_percentage: int = 50,
    output_folder: Optional[Union[str, os.PathLike]] = None,
) -> Optional[Path]:
    """Land-Sea mask operations.

    Parameters
    ----------
    ds: Union[xr.Dataset, str, os.PathLike]
    land_sea_mask: dict
    land_sea_percentage: int
    output_folder: str or os.PathLike, optional

    Returns
    -------
    Path
    """
    file_name = ""
    if isinstance(ds, (str, os.PathLike)):
        if output_folder is not None:
            output_folder = Path(output_folder)
            file_name = f"{Path(ds).stem}_land-sea-masked.nc"
        ds = xr.open_dataset(ds)

    if output_folder is not None and file_name == "":
        logging.warning(
            "Cannot generate filenames from xarray.Dataset objects. Consider writing NetCDF manually."
        )

    try:
        project = ds.attrs["project"]
    except KeyError:
        raise ValueError("No 'project' field found for given dataset.")

    if project in land_sea_mask.keys():
        logging.info(
            f"Masking variable with land-sea mask at {land_sea_percentage} % cutoff."
        )
        land_sea_mask_variable, lsm_file = land_sea_mask[project]
        lsm_raw = xr.open_dataset(lsm_file)
        try:
            lsm_raw = lsm_raw.rename({"longitude": "lon", "latitude": "lat"})
        except ValueError:
            raise

        lon_bounds = np.array([ds.lon.min(), ds.lon.max()])
        lat_bounds = np.array([ds.lat.min(), ds.lat.max()])

        lsm = subset.subset_bbox(
            lsm_raw,
            lon_bnds=lon_bounds,
            lat_bnds=lat_bounds,
        ).load()
        lsm = lsm.where(lsm[land_sea_mask_variable] > float(land_sea_percentage) / 100)
        if project == "era5":
            ds = ds.where(lsm[land_sea_mask].isel(time=0, drop=True).notnull())
            try:
                ds = ds.rename({"longitude": "lon", "latitude": "lat"})
            except ValueError:
                raise
        elif project in ["merra2", "cfsr"]:
            ds = ds.where(lsm[land_sea_mask].notnull())

        ds.attrs["land_sea_cutoff"] = f"{land_sea_percentage} %"

        if len(file_name) > 0:
            out = output_folder / file_name
            ds.to_netcdf(out)
            return out
        return ds
    raise RuntimeError(f"Project `{project}` was not found in land-sea masks.")


def delayed_write(
    ds: xr.Dataset,
    outfile: Path,
    target_chunks: dict,
    output_format: str,
    overwrite: bool,
):
    """

    Parameters
    ----------
    ds: Union[xr.Dataset, str, os.PathLike]
    outfile
    target_chunks
    output_format
    overwrite

    Returns
    -------

    """
    # Set correct chunks in encoding options
    kwargs = dict()
    kwargs["encoding"] = dict()
    for name, da in ds.data_vars.items():
        chunks = list()
        for dim in da.dims:
            if dim in target_chunks.keys():
                chunks.append(target_chunks[str(dim)])
            else:
                chunks.append(len(da[dim]))

        if output_format == "netcdf":
            kwargs["encoding"][name] = {
                "chunksizes": chunks,
                "zlib": True,
            }
            kwargs["compute"] = False
            if not overwrite:
                kwargs["mode"] = "a"
        elif output_format == "zarr":
            ds = ds.chunk(target_chunks)
            kwargs["encoding"][name] = {
                "chunks": chunks,
                "compressor": zarr.Blosc(),
            }
            kwargs["compute"] = False
            if overwrite:
                kwargs["mode"] = "w"
    if kwargs["encoding"]:
        kwargs["encoding"]["time"] = {"dtype": "int32"}

    return getattr(ds, f"to_{output_format}")(outfile, **kwargs)
