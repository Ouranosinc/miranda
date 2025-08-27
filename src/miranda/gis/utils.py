"""Utility functions for GIS operations."""

from __future__ import annotations
import datetime
import logging
import warnings

import numpy as np
import xarray as xr


logger = logging.getLogger("miranda.gis.utils")

__all__ = [
    "conservative_regrid",
    "threshold_mask",
]


def _simple_fix_dims(d: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
    """
    Adjust dimensions found in a file so that it can be used for regridding purposes.

    Parameters
    ----------
    d : xr.Dataset or xr.DataArray
        The dataset to adjust.

    Returns
    -------
    xr.Dataset or xr.DataArray
        The adjusted dataset.
    """
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


def conservative_regrid(ds: xr.DataArray | xr.Dataset, ref_grid: xr.DataArray | xr.Dataset) -> xr.DataArray | xr.Dataset:
    """
    Perform a conservative_normed regridding.

    Parameters
    ----------
    ds : xr.DataArray or xr.Dataset
        The dataset to regrid.
    ref_grid : xr.DataArray or xr.Dataset
        The reference grid.

    Returns
    -------
    xr.DataArray or xr.Dataset
        The regridded dataset.
    """
    try:
        import xesmf as xe  # noqa
    except ModuleNotFoundError:
        logger.error("This function requires the `xesmf` library which is not installed. Regridding step will be skipped.")
        raise

    ref_grid = _simple_fix_dims(ref_grid)
    method = "conservative_normed"

    msg = f"Performing regridding and masking with `xesmf` using method: {method}."
    logging.info(msg)

    regridder = xe.Regridder(ds, ref_grid, method, periodic=False)
    ds = regridder(ds)

    ds.attrs["history"] = f"{datetime.datetime.now()}:Regridded dataset using xesmf with method: {method}. {ds.attrs.get('history')}".strip()
    return ds


def threshold_mask(
    ds: xr.Dataset | xr.DataArray,
    *,
    mask: xr.Dataset | xr.DataArray,
    mask_cutoff: float | bool = False,
) -> xr.Dataset | xr.DataArray:
    """
    Land-Sea mask operations.

    Parameters
    ----------
    ds : xr.Dataset or str or os.PathLike
        The dataset to be masked.
    mask : xr.Dataset or xr.DataArray
        The land-sea mask.
    mask_cutoff : float or bool
        The mask cutoff value.

    Returns
    -------
    xr.Dataset or xr.DataArray
        The masked dataset.
    """
    mask = _simple_fix_dims(mask)

    if isinstance(mask, xr.Dataset):
        if len(mask.data_vars) == 1:
            mask_variable = list(mask.data_vars)[0]
            mask = mask[mask_variable]
        else:
            raise ValueError("More than one data variable found in land-sea mask. Supply a DataArray instead.")
    else:
        mask_variable = mask.name

    try:
        from clisops.core import subset_bbox  # noqa

        log_msg = f"Masking dataset with {mask_variable}."
        if mask_cutoff:
            log_msg = f"{log_msg.strip('.')} at `{mask_cutoff}` cutoff value."
        logging.info(log_msg)

        lon_bounds = np.array([ds.lon.min(), ds.lon.max()])
        lat_bounds = np.array([ds.lat.min(), ds.lat.max()])

        mask_subset = subset_bbox(
            mask,
            lon_bnds=lon_bounds,
            lat_bnds=lat_bounds,
        ).load()
    except ModuleNotFoundError:
        log_msg = "This function requires the `clisops` library which is not installed. subsetting step will be skipped."
        warnings.warn(log_msg, stacklevel=2)
        mask_subset = mask.load()

    if mask_subset.dtype == bool:
        if mask_cutoff:
            logging.warning("Mask value cutoff set for boolean mask. Ignoring.")
        mask_subset = mask_subset.where(mask)
    else:
        mask_subset = mask_subset.where(mask >= mask_cutoff)
    ds = ds.where(mask_subset.notnull())

    if mask_subset.min() >= 0:
        if mask_subset.max() <= 1.00000001:
            cutoff_info = f"{mask_cutoff * 100} %"
        elif mask_subset.max() <= 100.00000001:
            cutoff_info = f"{mask_cutoff} %"
        else:
            cutoff_info = f"{mask_cutoff}"
    else:
        cutoff_info = f"{mask_cutoff}"
    ds.attrs["mask_cutoff"] = cutoff_info

    prev_history = ds.attrs.get("history", "")
    history_msg = f"Mask calculated using `{mask_variable}`."
    if mask_cutoff:
        history_msg = f"{history_msg.strip('.')} with cutoff value `{cutoff_info}`."
    history = f"{history_msg} {prev_history}".strip()
    ds.attrs.update(dict(history=history))

    return ds
