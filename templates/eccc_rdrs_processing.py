import logging
from pathlib import Path

from miranda.convert.eccc_rdrs import convert_rdrs, rdrs_to_daily
from miranda.io import concat_rechunk_zarr


def main():
    home = Path("~").expanduser()
    dask_dir = home.joinpath("tmpout", "dask")
    dask_dir.mkdir(parents=True, exist_ok=True)
    dask_kwargs = dict(
        n_workers=8,
        threads_per_worker=4,
        memory_limit="7GB",
        dashboard_address=8998,
        local_directory=dask_dir,
        silence_logs=logging.ERROR,
    )

    # convert_rdrs(
    #     project="rdrs-v21",
    #     input_folder=Path(home).joinpath("RDRS_v21", "caspar"),
    #     output_folder=Path(home).joinpath("RDRS_v21", "tmp/ECCC/RDRS_v21/NAM"),
    #     output_format="zarr",
    #     working_folder=Path(home).joinpath("tmpout", "rdrs"),
    #     **dask_kwargs,
    # )

    rdrs_to_daily(
        project="rdrs-v21",
        input_folder=Path(home).joinpath(
            "RDRS_v21", "tmp/ECCC/RDRS_v21/NAM/1hr"
        ),
        output_folder=Path(home).joinpath("RDRS_v21", "tmp/ECCC/RDRS_v21/NAM/day"),
        working_folder=Path(home).joinpath("tmpout", "rdrs1"),
        overwrite=False,
        year_start=2002,
        **dask_kwargs,
    )

    # for freq in ["1hr", "day"]:
    #     infolder = Path(home).joinpath("RDRS_v2.1", f"tmp/ECCC/RDRS_v2.1/NAM/{freq}")
    #     for variable in [v for v in infolder.glob("*") if v.is_dir()]:
    #         concat_rechunk_zarr(
    #             input_folder=variable,
    #             output_folder=Path(home).joinpath(
    #                 "RDRS_v2.1", f"converted/ECCC/RDRS_v2.1/NAM/{freq}/{variable.name}"
    #             ),
    #             overwrite=False,
    #         )


if __name__ == "__main__":
    main()
