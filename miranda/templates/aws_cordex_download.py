import logging
from os import getenv
from pathlib import Path

import dask

from miranda.ncar import cordex_aws_download

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__file__)

dask.config.set(**{"array.slicing.split_large_chunks": True})

if __name__ == "__main__":

    out_files = getenv("out")

    target_folder = Path(out_files)

    search = dict(
        variable=["tasmax", "tasmin", "tas", "pr"],
        scenario=["rcp45", "rcp85", "hist"],
        grid="NAM-22i",
        bias_correction="raw",
    )

    cordex_aws_download(target_folder, search=search)
