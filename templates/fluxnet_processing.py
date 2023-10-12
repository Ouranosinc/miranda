import logging
from pathlib import Path

from miranda.convert.fluxnet import convert_fluxnet


def main():
    home = Path("~").expanduser()
    dask_dir = home.joinpath("tmpout", "dask")
    dask_dir.mkdir(parents=True, exist_ok=True)
    dask_kwargs = dict(
        n_workers=4,
        threads_per_worker=4,
        memory_limit="3GB",
        dashboard_address=8998,
        local_directory=dask_dir,
        silence_logs=logging.ERROR,
    )
    project = "fluxnet"
    convert_fluxnet(
        project=project,
        input_folder=Path(home).joinpath("fluxnet", "fluxnet2015_rawdata_subset"),
        output_folder=Path(home).joinpath("fluxnet", "FLUXNET2015_subset"),
        output_format="netcdf",
        working_folder=Path(home).joinpath("tmpout", "fluxnet"),
        overwrite=False,
        **dask_kwargs,
    )


if __name__ == "__main__":
    main()
