"""ECMWF TIGGE Conversion module."""

from __future__ import annotations

import itertools as it
import logging
import multiprocessing
import os
import shutil
import tempfile
from pathlib import Path

import xarray
from dask.diagnostics import ProgressBar

__all__ = ["tigge_convert"]


# FIXME: Is this function still pertinent?
def tigge_convert(
    source: os.PathLike | None = None,
    target: os.PathLike | None = None,
    processes: int = 8,
) -> None:
    """
    Convert TIGGE grib2 file to netCDF format.

    Parameters
    ----------
    source : os.PathLike, optional
        The source directory containing the TIGGE files.
    target : os.PathLike, optional
        The target directory to save the converted files.
    processes : int
        The number of processes to use for the conversion.
    """

    def _tigge_convert(fn):
        """
        Launch reformatting function.

        Parameters
        ----------
        fn : tuple
            The file and output folder.
        """
        infile, output_folder = fn
        try:
            for f in Path(infile.parent).glob(infile.name.replace(".grib", "*.idx")):
                f.unlink(missing_ok=True)

            ds = xarray.open_dataset(
                infile,
                engine="cfgrib",
                chunks="auto",
            )

            encoding = {var: dict(zlib=True) for var in ds.data_vars}
            encoding["time"] = {"dtype": "single"}
            tf = tempfile.NamedTemporaryFile(suffix=".nc", delete=False)

            with ProgressBar():
                msg = f"Converting `{infile.name}`."
                logging.info(msg)
                ds.to_netcdf(
                    tf.name, format="NETCDF4", engine="h5netcdf", encoding=encoding
                )

            shutil.move(
                tf.name,
                output_folder.joinpath(infile.name.replace(".grib", ".nc")).as_posix(),
            )

        except ValueError as err:
            msg = f"error converting {infile.name} : File may be corrupted : {err}"
            logging.error(msg)

    if source is None:
        source = Path().cwd().joinpath("downloaded")

    all_files = list(Path(source).glob("*.grib2"))
    if len(all_files) == 0:
        raise FileNotFoundError("TIGGE files not found.")

    if target is None:
        target = Path().cwd().joinpath("converted")
    else:
        target = Path(target)
    target.mkdir(exist_ok=True)

    with multiprocessing.Pool(processes=processes) as p:
        combs = list(it.product(*[all_files, [target]]))
        p.map(_tigge_convert, combs)
        p.close()
        p.join()
