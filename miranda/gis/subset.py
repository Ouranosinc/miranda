import logging
import logging.config
from pathlib import Path
from typing import List, Optional, Tuple, Union

import fiona
import geojson
import numpy as np
import rasterio.crs
import xarray

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["subset_shape", "subsetting_domains"]


def subsetting_domains(domain: str) -> np.array:
    """Provides the bounding box coordinates for specific domains.
    Parameters
    ----------
    domain: {"global", "nam", "can", "qc", "mtl"}
    Returns
    -------
    np.array
      North, West, South, and East coordinates
    """
    region = None

    if domain.upper() == "GLOBAL":
        region = np.array([90.0, -180.0, -90.0, 180.0])
    elif domain.upper() in ["AMNO", "NAM"]:
        region = np.array([90.0, -180.0, 10.0, -10.0])
    elif domain.upper() == "CAN":
        region = np.array([83.5, -141.0, 41.5, -52.5])
    elif domain.upper() == "QC":
        region = np.array([63.0, -80.0, 44.5, -57.0])
    elif domain.upper() == "MTL":
        region = np.array([45.75, -74.05, 45.3, -73.4])
    if region is not None:
        return region

    raise NotImplementedError(domain)


def _read_geometries(
    shape: Union[str, Path], crs: Optional[Union[str, int, dict]] = None
) -> Tuple[List[geojson.geometry.Geometry], rasterio.crs.CRS]:
    """
    A decorator to perform a check to verify a geometry is valid.
    Returns the function with geom set to the shapely Shape object.
    """
    try:
        if shape is None:
            raise ValueError
    except (KeyError, ValueError):
        logging.exception("No shape provided.")
        raise

    geom = list()
    geometry_types = list()
    try:
        with fiona.open(shape) as fio:
            logging.info("Vector read OK.")
            if crs:
                shape_crs = rasterio.crs.CRS.from_user_input(crs)
            else:
                shape_crs = rasterio.crs.CRS(fio.crs or 4326)
            for i, feat in enumerate(fio):
                g = geojson.GeoJSON(feat)
                geom.append(g["geometry"])
                geometry_types.append(g["geometry"]["type"])
    except fiona.errors.DriverError:
        logging.exception("Unable to read shape.")
        raise

    if len(geom):
        logging.info("Shapes found are {}.".format(", ".join(set(geometry_types))))
        return geom, shape_crs
    else:
        raise RuntimeError("No geometries found.")


def subset_shape(
    da: Union[xarray.DataArray, xarray.Dataset],
    shape: Union[str, Path],
    start_date: Optional[str] = None,
    end_date: Optional[str] = None,
    da_crs: Optional[str] = None,
    geometry: Optional[List[geojson.GeoJSON]] = None,
    shape_crs: Optional[str] = None,
) -> Union[xarray.DataArray, xarray.Dataset]:
    """Subset a DataArray or Dataset spatially (and temporally) using a vector shape and date selection.
    Return a subsetted data array for grid points falling within the area of a polygon and/or MultiPolygon shape,
      or grid points along the path of a LineString and/or MultiLineString.
    Parameters
    ----------
    da : Union[xarray.DataArray, xarray.Dataset]
      Input data.
    shape : Union[str, Path, geometry.GeometryCollection]
      Path to a single-layer vector file, or a shapely GeometryCollection object.
    start_date : Optional[str]
      Start date of the subset.
      Date string format -- can be year ("%Y"), year-month ("%Y-%m") or year-month-day("%Y-%m-%d").
      Defaults to first day of input data-array.
    end_date : Optional[str]
      End date of the subset.
      Date string format -- can be year ("%Y"), year-month ("%Y-%m") or year-month-day("%Y-%m-%d").
      Defaults to last day of input data-array.
    geometry: Optional[List[geojson.GeoJSON]]
      A list of all GeoJSON shapes to be used in clipping the xarray.DataArray or xarray.Dataset
    da_crs : Optional[Union[int, dict, str]]
      CRS of the xarray.DataArray or xarray.Dataset provided. Default: dict(epsg=4326).
    shape_crs : Optional[Union[int, dict, str]]
      CRS of the geometries provided. If passing GeometryCollections as shapes, CRS must be explicitly stated.
    Returns
    -------
    Union[xarray.DataArray, xarray.Dataset]
      Subsetted xarray.DataArray or xarray.Dataset
    Warnings
    --------
    This functions relies on the rioxarray library and requires xarray Datasets and DataArrays that have been read with
     the `rioxarray` library imported. Attempting to use this function with pure xarray objects will raise exceptions.
    Examples
    --------
    >>> from miranda import gis
    >>> import rioxarray
    >>> import xarray as xr
    >>> ds = xarray.open_dataset('pr.day.nc')
    Subset lat lon and years
    >>> prSub = gis.subset_shape(ds.pr, shape="/path/to/polygon.shp", start_yr='1990', end_yr='1999')
    Subset data array lat, lon and single year
    >>> prSub = gis.subset_shape(ds.pr, shape="/path/to/polygon.shp", start_yr='1990', end_yr='1990')
    Subset data array single year keep entire lon, lat grid
    >>> prSub = gis.subset_bbox(ds.pr, start_yr='1990', end_yr='1990') # one year only entire grid
    Subset multiple variables in a single dataset
    >>> ds = xarray.open_mfdataset(['pr.day.nc','tas.day.nc'])
    >>> dsSub = gis.subset_bbox(ds, shape="/path/to/polygon.shp", start_yr='1990', end_yr='1999')
     # Subset with year-month precision - Example subset 1990-03-01 to 1999-08-31 inclusively
    >>> prSub = gis.subset_time(ds.pr, shape="/path/to/polygon.shp", start_date='1990-03', end_date='1999-08')
    # Subset with specific start_dates and end_dates
    >>> prSub = \
            subset.subset_time(ds.pr, shape="/path/to/polygon.shp", start_date='1990-03-13', end_date='1990-08-17')
    """

    if geometry and shape_crs:
        shape_crs = rasterio.crs.CRS.from_user_input(shape_crs)
    else:
        geometry, shape_crs = _read_geometries(shape, crs=shape_crs)

    try:
        # NetCDF data doesn't typically have defined CRS. Ensure this is the case and append one if needed.
        if da.rio.crs is None:
            if da_crs is None:

                if "rlon" in da.dims or "rlat" in da.dims:
                    raise NotImplementedError("Rotated poles are not supported.")

                else:
                    crs = rasterio.crs.CRS.from_epsg(4326)
                    if np.any(da.lon < -180) or np.any(da.lon > 360):
                        raise rasterio.crs.CRSError(
                            "NetCDF doesn't seem to be in EPSG:4326. Set CRS manually."
                        )

                    # Convert longitudes from 0,+360 to -180,+180
                    if np.any(da.lon > 180):
                        lon_attrs = da.lon.attrs.copy()
                        fix_lon = da.lon.values
                        fix_lon[fix_lon > 180] = fix_lon[fix_lon > 180] - 360

                        # Correct lon_bnds in xarray.Datasets
                        if isinstance(da, xarray.Dataset):
                            if "lon_bnds" in da.data_vars:
                                fix_lon_bnds = da.lon_bnds.values
                                fix_lon_bnds[fix_lon_bnds > 180] = (
                                    fix_lon_bnds[fix_lon_bnds > 180] - 360
                                )
                        da = da.assign_coords(lon=fix_lon)
                        da = da.sortby("lon")
                        da.lon.attrs = lon_attrs
            else:
                crs = rasterio.crs.CRS.from_user_input(da_crs)
        else:
            crs = da.rio.crs

    except Exception as e:
        logging.exception(e)
        raise

    if shape_crs != crs:
        raise rasterio.crs.CRSError("Shape and Raster CRS are not the same.")

    if isinstance(da, xarray.Dataset):
        # Create a new empty xarray.Dataset and populate with corrected/clipped variables
        ds_out = xarray.Dataset(data_vars=None, attrs=da.attrs)
        for v in da.data_vars:
            if "lon" in da[v].dims and "lat" in da[v].dims:
                dss = da[v]
                # Identify spatial dimensions
                dss.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
                dss.rio.write_crs(crs, inplace=True)
                ds_out[v] = dss.rio.clip(
                    geometry, crs=dss.rio.crs, all_touched=True, drop=True, invert=False
                )

        for v in da.data_vars:
            if not ("lon" in da[v].dims and "lat" in da[v].dims):
                if "lat" in da[v].dims:
                    ds_out[v] = da[v].sel(lat=ds_out.lat)
                elif "lon" in da[v].dims:
                    ds_out[v] = da[v].sel(lon=ds_out.lon)
                else:
                    ds_out[v] = da[v]
    else:
        da.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        da.rio.write_crs(crs, inplace=True)
        ds_out = da.rio.clip(
            geometry, crs=crs, all_touched=True, drop=True, invert=False
        )

    return ds_out
