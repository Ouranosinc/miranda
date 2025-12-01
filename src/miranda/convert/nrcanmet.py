"""NRCANmet (ANUSPLIN) interpolated station data conversion module."""
from __future__ import annotations
from pathlib import Path

import xarray as xr

from miranda.convert._data_corrections import (
    dataset_conversion,
)


__all__ = ["convert_nrcanmet"]


def convert_nrcanmet(infile: str | Path, engine: str = "h5netcdf") -> xr.Dataset:
    """
    Convert the NRCanMET netCDF files to production-ready CF-compliant netCDFs.

    Parameters
    ----------
    infolder : str or Path
        The path to the NRCanMET netCDF files.
    outfolder : str or Path
        The output directory.
    """
    if isinstance(infile, str):
        infile = Path(infile)

    ds = xr.open_dataset(infile, chunks={}, engine=engine)
    if "crs" in ds.data_vars:
        ds = ds.assign_coords(crs=ds.crs)
    ds = dataset_conversion(ds, project="NRCanMET")
    return ds
