import logging
import os
from pathlib import Path

from miranda.convert.eccc_rdrs import concat_zarr, convert_rdrs, rdrs_to_daily

home = os.path.expanduser("~")
dask_dir = Path(home).joinpath("tmpout", "dask")
dask_dir.mkdir(parents=True, exist_ok=True)
dask_kwargs = dict(
    n_workers=5,
    threads_per_worker=5,
    memory_limit="7GB",
    dashboard_address=8999,
    local_directory=dask_dir,
    silence_logs=logging.ERROR,
)

project = "rdrs-v2.1"
convert_rdrs(
    project=project,
    input_folder=Path(home).joinpath("RDRS_v2.1", "caspar"),
    output_folder=Path(home).joinpath("RDRS_v2.1", "tmp/ECCC/RDRS_v2.1/NAM"),
    output_format="zarr",
    working_folder=Path(home).joinpath("tmpout", "rdrs"),
    **dask_kwargs,
)

rdrs_to_daily(
    input_folder=Path(home).joinpath("RDRS_v2.1", "converted/ECCC/RDRS_v2.1/NAM/1hr"),
    output_folder=Path(home).joinpath("RDRS_v2.1", "tmp/ECCC/RDRS_v2.1/NAM/day"),
    working_folder=Path(home).joinpath("tmpout", "rdrs"),
    overwrite=False,
)


for freq in ["1hr", "day"]:
    infolder = Path(home).joinpath("RDRS_v2.1", f"tmp/ECCC/RDRS_v2.1/NAM/{freq}")
    for vv in [i for i in infolder.glob("*") if i.is_dir()]:
        concat_zarr(
            input_folder=vv,
            output_folder=Path(home).joinpath(
                "RDRS_v2.1", f"converted/ECCC/RDRS_v2.1/NAM/{freq}/{vv.name}"
            ),
            overwrite=False,
        )
