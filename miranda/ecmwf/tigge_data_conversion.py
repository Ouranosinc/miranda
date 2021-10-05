import itertools as it
import logging
import multiprocessing
import os
import shutil
import tempfile
from datetime import datetime as dt
from pathlib import Path

import xarray
from dask.diagnostics import ProgressBar

# logging file configuration
logging.basicConfig(
    filename="{}_{}.log".format(dt.now().strftime("%Y%m%d"), Path(__file__).stem),
    level=logging.INFO,
    datefmt="%H:%M:%S",
)

# set up logging to console
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# format output to the console
formatter = logging.Formatter("%(name)s : %(asctime)s :  %(levelname)s : %(message)s")
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger("").addHandler(console)
logger = logging.getLogger(__name__)


def main(source: os.PathLike, target: os.PathLike):

    all_files = Path(source).glob("*.grib2")

    target = Path(target)
    target.mkdir(exist_ok=True)

    p = multiprocessing.Pool(processes=10)

    combs = list(it.product(*[all_files, [target]]))
    p.map(convert_tigge, combs)
    p.close()
    p.join()


def convert_tigge(fn):
    """Convert grib2 file to netCDF format."""
    infile, outfolder = fn
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
            outfolder.joinpath(infile.name.replace(".grib", ".nc")).as_posix(),
        )

    except ValueError:
        print(f"error converting {infile.name} : File may be corrupted")
        pass


if __name__ == "__main__":
    source = Path().cwd().joinpath("download")
    target = Path().cwd().joinpath("converted")
    main(source, target)
