from pathlib import Path
from tempfile import tempdir

from dask.distributed import Client

from miranda.ecmwf import rechunk_ecmwf

if __name__ == "__main__":
    step = "hourly"  # "daily
    target_project = "era5-land"  # "era5-single-levels"
    outfmt = "nc"  # "zarr"

    base_path = Path(f"/path/to/{target_project}/downloaded/")
    new_path = Path(f"/path/to/{target_project}/rechunked/")

    with Client(
        n_workers=2,
        threads_per_worker=1,
        dashboard_address=8786,
        memory_limit="4GB",
        local_directory=Path(tempdir),
    ):
        rechunk_ecmwf(
            project=target_project,
            input_folder=base_path,
            output_folder=new_path,
            time_step=step,
            output_format=outfmt,
        )
