from __future__ import annotations

import logging.config

import numpy as np
import xarray as xr

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["subset_domain", "subsetting_domains", "add_ar6_regions"]

_gis_import_error_message = (
    "`{}` requires installation of the miranda GIS libraries. These can be installed using the"
    " `pip install miranda[gis]` recipe or via Anaconda (`conda env update -n miranda-env -f environment.yml`)"
    " from the miranda repository source files."
)


def subset_domain(
    ds: xr.Dataset | xr.DataArray, domain: str, **kwargs
) -> xr.Dataset | xr.DataArray:
    r"""Subset an xarray object according to a specific domain.

    Notes
    -----
    Requires installation of GIS libraries.

    Parameters
    ----------
    ds: xarray.Dataset or xarray.DataArray
    domain: str
    \*\*kwargs

    Returns
    -------
    xarray.Dataset or xarray.DataArray
    """
    try:
        from clisops.core.subset import subset_bbox  # noqa
    except ModuleNotFoundError:
        msg = _gis_import_error_message.format(subset_domain.__name__)
        raise ModuleNotFoundError(msg)

    region = subsetting_domains(domain)
    lon_values = np.array([region[1], region[3]])
    lat_values = np.array([region[0], region[2]])

    ds = subset_bbox(ds, lon_bnds=lon_values, lat_bnds=lat_values, **kwargs)

    return ds


def subsetting_domains(domain: str) -> list:
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
# def _read_geometries(
#     shape: Union[str, Path], crs: Optional[Union[str, int, dict]] = None
# ) -> Tuple[List[geojson.geometry.Geometry], Any]:
#     """Perform a check to verify a geometry is valid.
#
#     Returns the function with geom set to the shapely Shape object.
#     """
#     try:
#         from pyproj import CRS
#     except ModuleNotFoundError:
#         msg = gis_import_error_message.format(add_ar6_regions.__name__)
#         raise ModuleNotFoundError(msg)
#
#     try:
#         if shape is None:
#             raise ValueError
#     except (KeyError, ValueError):
#         logging.exception("No shape provided.")
#         raise
#
#     geom = list()
#     geometry_types = list()
#     try:
#         with fiona.open(shape) as fio:
#             logging.info("Vector read OK.")
#             if crs:
#                 shape_crs = CRS.from_user_input(crs)
#             else:
#                 shape_crs = CRS(fio.crs or 4326)
#             for i, feat in enumerate(fio):
#                 g = geojson.GeoJSON(feat)
#                 geom.append(g["geometry"])
#                 geometry_types.append(g["geometry"]["type"])
#     except fiona.errors.DriverError:
#         logging.exception("Unable to read shape.")
#         raise
#
#     if len(geom) > 0:
#         logging.info(f"Shapes found are: {', '.join(set(geometry_types))}.")
#         return geom, shape_crs
#     raise RuntimeError("No geometries found.")


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
        msg = _gis_import_error_message.format(add_ar6_regions.__name__)
        raise ImportError(msg)

    mask = regionmask.defined_regions.ar6.all.mask(ds.lon, ds.lat)
    ds = ds.assign_coords(region=mask)
    return ds
