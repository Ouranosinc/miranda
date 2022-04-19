from os import getenv
from pathlib import Path
from tempfile import tempdir

from dask.distributed import Client

from miranda.convert import rechunk_reanalysis

if __name__ == "__main__":
    step = "hourly"  # "daily
    target_project = "era5-land"  # "era5-single-levels"
    outfmt = "nc"  # "zarr"

    in_files = getenv("in")
    out_files = getenv("out")

    input_path = Path(in_files)
    output_path = Path(out_files)

    with Client(
        n_workers=2,
        threads_per_worker=1,
        dashboard_address=8786,
        memory_limit="4GB",
        local_directory=Path(tempdir),
    ):
        rechunk_reanalysis(
            project=target_project,
            input_folder=input_path,
            output_folder=output_path,
            time_step=step,
            output_format=outfmt,
        )
