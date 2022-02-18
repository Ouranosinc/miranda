import logging
import logging.config
from pathlib import Path
from typing import List, Optional, Tuple, Union

import fiona
import geojson
import numpy as np
import rasterio.crs

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["subsetting_domains"]


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

    if len(geom) > 0:
        logging.info(f"Shapes found are: {', '.join(set(geometry_types))}.")
        return geom, shape_crs
    else:
        raise RuntimeError("No geometries found.")
