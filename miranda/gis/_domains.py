import logging
import logging.config
from pathlib import Path
from typing import List, Optional, Tuple, Union

import fiona
import geojson
import numpy as np
import regionmask
import xarray as xr
from clisops.core.subset import subset_bbox
from pyproj.crs import CRS

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["subset_domain", "subsetting_domains", "add_ar6_regions"]


def subset_domain(
    ds: Union[xr.Dataset, xr.DataArray], domain: str, **kwargs
) -> Union[xr.Dataset, xr.DataArray]:
    region = subsetting_domains(domain)
    lon_values = np.array([region[1], region[3]])
    lat_values = np.array([region[0], region[2]])

    ds = subset_bbox(ds, lon_bnds=lon_values, lat_bnds=lat_values, **kwargs)

    return ds


def subsetting_domains(domain: str) -> List:
    """Provides the bounding box coordinates for specific domains.

    Parameters
    ----------
    domain : {"global", "nam", "can", "qc", "mtl"}

    Returns
    -------
    np.array
      North, West, South, and East coordinates
    """
    region = None

    if domain.upper() == "GLOBAL":
        region = [90.0, -180.0, -90.0, 180.0]
    elif domain.upper() in ["AMNO", "NAM"]:
        region = [90.0, -179.9, 10.0, -10.0]
    elif domain.upper() == "CAN":
        region = [83.5, -141.0, 41.5, -52.5]
    elif domain.upper() == "QC":
        region = [63.0, -80.0, 44.5, -57.0]
    elif domain.upper() == "MTL":
        region = [45.75, -74.05, 45.3, -73.4]
    if region is not None:
        return region

    raise NotImplementedError(domain)


# FIXME: This approach will not work with Fiona 1.9+
def _read_geometries(
    shape: Union[str, Path], crs: Optional[Union[str, int, dict]] = None
) -> Tuple[List[geojson.geometry.Geometry], CRS]:
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
                shape_crs = CRS.from_user_input(crs)
            else:
                shape_crs = CRS(fio.crs or 4326)
            for i, feat in enumerate(fio):
                g = geojson.GeoJSON(feat)
                geom.append(g["geometry"])
                geometry_types.append(g["geometry"]["type"])
    except fiona.errors.DriverError:
        logging.exception("Unable to read shape.")
        raise

    if len(geom) > 0:
        logging.info(f"Shapes found are: {', '.join(set(geometry_types))}.")
        return geom, shape_crs
    raise RuntimeError("No geometries found.")


def add_ar6_regions(ds: xr.Dataset) -> xr.Dataset:
    """Add the IPCC AR6 Regions to dataset.

    Parameters
    ----------
    ds : xarray.Dataset

    Returns
    -------
    xarray.Dataset
    """

    mask = regionmask.defined_regions.ar6.all.mask(ds.lon, ds.lat)
    ds = ds.assign_coords(region=mask)
    return ds
