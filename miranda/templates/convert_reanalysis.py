from os import getenv
from pathlib import Path

import dask

if __name__ == "__main__":
    dask.config.set(local_directory=f"{Path(__file__).parent}/dask_workers/")

    in_files = getenv("in")
    out_files = getenv("out")

    input_path = Path(in_files)
    output_path = Path(out_files)

    variables = ["tas", "pr"]
    start, end = "1950", "2021"
    domains = None  # {"QC", "CAN", "AMNO"}
