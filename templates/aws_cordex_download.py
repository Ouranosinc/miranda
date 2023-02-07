import logging
from os import getenv
from pathlib import Path

import dask

from miranda.ncar._aws_cordex import cordex_aws_download  # noqa

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__file__)

dask.config.set(**{"array.slicing.split_large_chunks": True})

if __name__ == "__main__":
    out_files = getenv("out")
    target_folder = Path(out_files).expanduser()
    domain = None  # "AMNO"

    search = dict(
        variable=["tasmax", "tasmin", "tas", "pr"],
        scenario=["rcp45", "rcp85", "hist"],
        grid=["NAM-22i"],
        bias_correction=["raw"],
    )

    cordex_aws_download(target_folder, search=search, domain=domain)
