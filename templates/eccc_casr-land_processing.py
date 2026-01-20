import logging
from pathlib import Path

from miranda.convert.eccc_rdrs import convert_rdrs, rdrs_to_daily
from miranda.io import concat_rechunk_zarr


def main():
    vars_to_process = None #["tas", "pr"]
    project = "casr-land-v21"

    home = Path("~").expanduser()
    dask_dir = home.joinpath("tmpout", project, "dask")
    dask_dir.mkdir(parents=True, exist_ok=True)
    dask_kwargs = dict(
        n_workers=6,
        threads_per_worker=6,
        memory_limit="12GB",
        dashboard_address=8991,
        local_directory=dask_dir,
        silence_logs=logging.ERROR,
    )


    convert_rdrs(
        project=project,
        input_folder=Path(home).joinpath("data", "casr-land", "raw", "final"),
        output_folder=Path(home).joinpath("data", "casr-land", "processed"),
        output_format="zarr",
        working_folder=Path(home).joinpath("tmpout", project, "working"),
        cfvariable_list=vars_to_process,
        year_start=1981,
        year_end=1981,
        overwrite=False,
        **dask_kwargs,
    )

    rdrs_to_daily(
        project=project,
        input_folder=Path(home).joinpath("data", "casr-land", "processed", "1hr"),
        output_folder=Path(home).joinpath("data", "casr-land", "processed", "day"),
        working_folder=Path(home).joinpath("tmpout", project, "working"),
        overwrite=False,
        year_start=None,
        year_end=None,
        process_variables=vars_to_process,
        **dask_kwargs,
    )
    # #

    # for freq in ["day", "1hr"]:
    #     infolder = Path(home).joinpath("RDRS_v21", f"tmp/ECCC/RDRS_v21/NAM/{freq}")
    #     for variable in [
    #         v for v in infolder.glob("*") if v.is_dir() and v.name in vars_to_process
    #     ]:
    #         concat_rechunk_zarr(
    #             project=project,
    #             freq=freq,
    #             input_folder=variable,
    #             output_folder=Path(home).joinpath(
    #                 "RDRS_v21", f"converted/ECCC/RDRS_v21/NAM/{freq}/{variable.name}"
    #             ),
    #             overwrite=False,
    #             **dask_kwargs,
    #         )


if __name__ == "__main__":
    main()
