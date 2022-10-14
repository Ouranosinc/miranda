import logging.config
import os
from pathlib import Path
from typing import Dict, Optional, Union

import netCDF4
import numpy as np
import xarray as xr
import zarr
from clisops.core import subset
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "get_chunks_on_disk",
    "add_ar6_regions",
    "daily_aggregation",
    "delayed_write",
]


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
