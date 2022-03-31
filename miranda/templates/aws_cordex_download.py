import logging
from pathlib import Path

import dask

from miranda.ncar import cordex_aws_download

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__file__)

dask.config.set(**{"array.slicing.split_large_chunks": True})

if __name__ == "__main__":
    search = dict(
        variable=["tasmax", "tasmin", "tas", "pr"],
        scenario=["rcp45", "rcp85", "hist"],
        grid="NAM-22i",
        bias_correction="raw",
    )
    target_folder = Path("/tank/smith/cordex_ncar_interpolated")

    cordex_aws_download(target_folder, search=search)
