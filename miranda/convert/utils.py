import hashlib
import logging.config
import os
import re
from pathlib import Path
from typing import Callable, Dict, Optional, Union

import netCDF4
import numpy as np
import xarray as xr
import zarr
from clisops.core import subset
from dask.delayed import delayed
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "add_ar6_regions",
    "daily_aggregation",
    "delayed_write",
    "find_version_hash",
    "get_chunks_on_disk",
]


def get_chunks_on_disk(file: Union[os.PathLike, str]) -> dict:
    """

    Parameters
    ----------
    file : str or os.PathLike
        File to be examined. Supports NetCDF and Zarr.

    Returns
    -------
    dict
    """
    chunks = dict()
    file = Path(file)

    if file.suffix.lower() in [".nc", ".nc4"]:
        with netCDF4.Dataset(file) as ds:
            for v in ds.variables:
                chunks[v] = dict()
                for ii, dim in enumerate(ds[v].dimensions):
                    chunks[v][dim] = ds[v].chunking()[ii]
    elif file.suffix.lower() == "zarr" and file.is_dir():
        with zarr.open(file, "r") as ds:  # noqa
            for v in ds.arrays():
                # Check if variable is chunked
                if v[1]:
                    chunks[v[0]] = v[1]
    else:
        raise NotImplementedError(f"File type: {file.suffix}.")
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


def daily_aggregation(ds: xr.Dataset) -> Dict[str, xr.Dataset]:
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
            "hur",
            "hus",
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
    ds: Union[xr.Dataset, xr.DataArray],
    *,
    land_sea_mask: xr.DataArray,
    land_sea_cutoff: float = 0.5,
) -> Union[xr.Dataset, xr.DataArray]:
    """Land-Sea mask operations.

    Parameters
    ----------
    ds : Union[xr.Dataset, str, os.PathLike]
    land_sea_mask : Union[xr.Dataset, xr.DataArray]
    land_sea_cutoff : int

    Returns
    -------
    Union[xr.Dataset, xr.DataArray]
    """
    logging.info(
        f"Masking variable with land-sea mask at `{land_sea_cutoff}` cutoff value."
    )
    if "lon" not in land_sea_mask.dims or "lat" not in land_sea_mask.dims:
        land_sea_mask = land_sea_mask.rename(
            {
                "longitude": "lon",
                "latitude": "lat",
                "Longitude": "lon",
                "Latitude": "lat",
                "lons": "lon",
                "lats": "lat",
            }
        )
    if "time" in land_sea_mask.dims:
        land_sea_mask = land_sea_mask.isel(time=0, drop=True)

    if isinstance(land_sea_mask, xr.Dataset):
        if len(land_sea_mask.data_vars) == 1:
            land_sea_mask = ds[list(ds.data_vars)[0]]
        else:
            raise ValueError(
                "More than one data variable found in land-sea mask. Supply a DataArray instead."
            )

    lon_bounds = np.array([ds.lon.min(), ds.lon.max()])
    lat_bounds = np.array([ds.lat.min(), ds.lat.max()])

    lsm = subset.subset_bbox(
        land_sea_mask,
        lon_bnds=lon_bounds,
        lat_bnds=lat_bounds,
    ).load()

    lsm = lsm.where(land_sea_mask > land_sea_cutoff)
    ds = ds.where(lsm.notnull())

    if lsm.min() >= 0 and lsm.max() <= 1:
        ds.attrs["land_sea_cutoff"] = f"{land_sea_cutoff * 100} %"
    elif lsm.min() >= 0 and lsm.max() <= 100:
        ds.attrs["land_sea_cutoff"] = f"{land_sea_cutoff} %"
    else:
        ds.attrs["land_sea_cutoff"] = f"{land_sea_cutoff}"

    return ds


def find_version_hash(file: Union[os.PathLike, str]) -> Dict:
    def _get_hash(f):
        hash_sha256_writer = hashlib.sha256()
        with open(f, "rb") as f_opened:
            hash_sha256_writer.update(f_opened.read())
        sha256sum = hash_sha256_writer.hexdigest()
        logging.info(f"Calculated sha256sum (starting: {sha256sum[:6]})")
        del hash_sha256_writer
        return sha256sum

    version_info = dict()
    possible_version = Path(file).parent.name
    if re.match(r"^v\d+", possible_version, re.IGNORECASE):
        version_info["version"] = Path(file).parent.name
        version_info["sha256sum"] = _get_hash(file)

    else:
        file_identity = str(Path(file).name).split(".")[0]
        possible_version_signature = Path(file).parent.glob(f"{file_identity}.*")
        for sig in possible_version_signature:
            found_version = re.search(r"\.(v\d+.+)$", sig.name, re.IGNORECASE)
            if found_version:
                try:
                    version_info["version"] = found_version.group()
                    version_info["sha256sum"] = int(sig.open().read())
                except ValueError:
                    continue
                break
        else:
            version_info["version"] = "vNotFound"
            version_info["sha256sum"] = _get_hash(file)

    return version_info


def delayed_write(
    ds: xr.Dataset,
    outfile: Union[str, os.PathLike],
    output_format: str,
    overwrite: bool,
    target_chunks: Optional[dict] = None,
) -> delayed:
    """

    Parameters
    ----------
    ds : Union[xr.Dataset, str, os.PathLike]
    outfile : str or os.PathLike
    target_chunks : dict
    output_format : {"netcdf", "zarr"}
    overwrite : bool

    Returns
    -------
    dask.delayed.delayed
    """
    # Set correct chunks in encoding options
    kwargs = dict()
    kwargs["encoding"] = dict()
    for name, da in ds.data_vars.items():
        chunks = list()
        for dim in da.dims:
            if target_chunks:
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
