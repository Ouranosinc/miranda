import itertools as it
import multiprocessing
import os
import shutil
import tempfile
from pathlib import Path
from typing import Optional

import xarray
from dask.diagnostics import ProgressBar


def tigge_convert(
    source: Optional[os.PathLike] = None,
    target: Optional[os.PathLike] = None,
    processes: int = 8,
) -> None:
    """Convert grib2 file to netCDF format.

    Parameters
    ----------
    source : os.PathLike, optional
    target : os.PathLike, optional
    processes : int

    Returns
    -------
    None
    """

    def _tigge_convert(fn):
        """Launch reformatting function."""
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
            tmpfile = tempfile.NamedTemporaryFile(suffix=".nc", delete=False)

            with ProgressBar():
                print("converting ", infile.name)
                ds.to_netcdf(
                    tmpfile.name, format="NETCDF4", engine="netcdf4", encoding=encoding
                )

            shutil.move(
                tmpfile.name,
                output_folder.joinpath(infile.name.replace(".grib", ".nc")).as_posix(),
            )

        except ValueError:
            print(f"error converting {infile.name} : File may be corrupted")

    if source is None:
        source = Path().cwd().joinpath("download")
    if target is None:
        target = Path().cwd().joinpath("converted")

    all_files = Path(source).glob("*.grib2")

    target = Path(target)
    target.mkdir(exist_ok=True)

    p = multiprocessing.Pool(processes=processes)

    combs = list(it.product(*[all_files, [target]]))
    p.map(_tigge_convert, combs)
    p.close()
    p.join()
