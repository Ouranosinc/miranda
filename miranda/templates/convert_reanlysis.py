from pathlib import Path

import dask

dask.config.set(local_directory=f"{Path(__file__).parent}/dask_workers/")


outfiles = Path("/path/to/outfiles")
variables = ["tas", "pr"]
start, end = "1950", "2021"
domains = None  # {"QC", "CAN", "AMNO"}
