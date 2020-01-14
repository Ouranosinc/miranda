#!/bin/env python3
import datetime
import logging
from logging import config
from pathlib import Path

import geopandas as gpd
import xarray as xr

from miranda.scripting import LOGGING_CONFIG

config.dictConfig(LOGGING_CONFIG)

try:
    import dask
    from dask.diagnostics import ProgressBar
    from dask.distributed import Client

    from multiprocessing.pool import ThreadPool

    dask.config.set(scheduler="threads", pool=ThreadPool(10))

except ImportError as e:
    dask = None
    ProgressBar = None
    ThreadPool = None
    msg = "{}: Modules not found. Continuing with Xarray only.".format(e)
    logging.warning(msg)

__all__ = ["subset_bbox", "subset_bbox_years", "subset_mask"]


def subset_bbox_years(
    nc_file,
    lat_bounds,
    lon_bounds,
    start_year: (int, bool) = None,
    end_year: (int, bool) = None,
    output=None,
):
    if output is None:
        output = Path(nc_file).parent
    elif not Path(output).exists():
        Path.mkdir(output, parents=True)

    outfile = Path.joinpath(output, Path(nc_file).basename.replace(".nc", "_subset.nc"))
    ds = xr.open_dataset(nc_file, chunks={"lat": 100, "lon": 100, "time": 10})

    ds_sub = ds.where(
        (ds.lon >= min(lon_bounds))
        & (ds.lon <= max(lon_bounds))
        & (ds.lat >= min(lat_bounds))
        & (ds.lat <= max(lat_bounds)),
        drop=True,
    )

    if start_year is not None & end_year is not None:
        ds_sub = ds.where(
            (ds.time.dt.year > start_year) & (ds.time.dt.year <= end_year), drop=True
        )

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds_sub.to_netcdf(outfile, format="NETCDF4", encoding=encoding)

    now = datetime.datetime.now()
    logging.info([outfile, now.strftime("%Y-%m-%d %X")])

    try:
        with ProgressBar():
            ds_sub.compute()
    except ModuleNotFoundError:
        ds_sub.compute()
    return outfile


def subset_bbox(ncfile, lat_bnds, lon_bnds, output):
    if not Path(output).exists():
        Path.mkdir(output, parents=True)

    outfile = Path.joinpath(output, Path(ncfile).basename.replace(".nc", "_subset.nc"))

    ds = xr.open_dataset(ncfile).chunk({"time": 10})

    ds_sub = ds.where(
        (ds.lon >= min(lon_bnds))
        & (ds.lon <= max(lon_bnds))
        & (ds.lat >= min(lat_bnds))
        & (ds.lat <= max(lat_bnds)),
        drop=True,
    )
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds_sub.to_netcdf(outfile, format="NETCDF4", encoding=encoding, compute=False)

    now = datetime.datetime.now()
    logging.info([outfile, now.strftime("%Y-%m-%d %X")])

    with ProgressBar():
        ds_sub.compute()


def subset_mask(dataset: xr.Dataset, polygon: gpd.GeoDataFrame) -> xr.DataArray:
    import numpy as np
    import pandas as pd
    from shapely.geometry import Point

    dataset = dataset.drop(dataset.data_vars)
    dataset = dataset.drop("time")
    if len(dataset.lon.shape) == 1 & len(dataset.lat.shape) == 1:
        # Create a 2d grid of lon, lat values
        lon1, lat1 = np.meshgrid(
            np.asarray(dataset.lon.values), np.asarray(dataset.lat.values)
        )
    else:
        lon1 = dataset.lon.values
        lat1 = dataset.lat.values

    # Create Pandas DataFrame from NetCDF lat lon points
    df = pd.DataFrame(
        {"id": np.arange(0, lon1.size), "lon": lon1.flatten(), "lat": lat1.flatten()}
    )

    df["Coordinates"] = list(zip(df.lon, df.lat))
    df["Coordinates"] = df["Coordinates"].apply(Point)

    # Create GeoDataFrame (spatially referenced)
    gdf_pts = gpd.GeoDataFrame(df, geometry="Coordinates")
    gdf_pts.crs = {"init": "epsg:4326"}

    # Spatial join GeoData points with region polygons
    point_in_poly = gpd.tools.sjoin(gdf_pts, polygon, how="left", op="within")

    # Extract polygon ids for points
    mask = point_in_poly["index_right"]

    polygon_mask = np.array(mask).reshape(lat1.shape[0], lat1.shape[1])
    polygon_mask = xr.DataArray(polygon_mask, coords=dataset.coords, dims=dataset.dims)

    return polygon_mask
