import logging
from pathlib import Path

from miranda.convert.eccc_rdrs import convert_rdrs, rdrs_to_daily
from miranda.io import concat_rechunk_zarr


def main():
    home = Path("~").expanduser()
    dask_dir = home.joinpath("tmpout", "dask")
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
        input_folder=Path(home).joinpath(
            "RDRS_v2.1", "converted/ECCC/RDRS_v2.1/NAM/1hr"
        ),
        output_folder=Path(home).joinpath("RDRS_v2.1", "tmp/ECCC/RDRS_v2.1/NAM/day"),
        working_folder=Path(home).joinpath("tmpout", "rdrs"),
        overwrite=False,
        **dask_kwargs,
    )

    for freq in ["1hr", "day"]:
        infolder = Path(home).joinpath("RDRS_v2.1", f"tmp/ECCC/RDRS_v2.1/NAM/{freq}")
        for variable in [v for v in infolder.glob("*") if v.is_dir()]:
            concat_zarr(
                input_folder=variable,
                output_folder=Path(home).joinpath(
                    "RDRS_v2.1", f"converted/ECCC/RDRS_v2.1/NAM/{freq}/{variable.name}"
                ),
                overwrite=False,
            )


if __name__ == "__main__":
    main()
