from os import getenv
from pathlib import Path

from miranda.eccc import (
    aggregate_stations,
    convert_flat_files,
    merge_converted_variables,
)

if __name__ == "__main__":
    time_step = "daily"
    n_workers = 3
    var_codes = [
        1,
        2,
        3,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    ]

    in_files = getenv("in")
    source_data = Path(in_files)

    origin_files = source_data.joinpath("source")

    daily = source_data.joinpath("daily")
    output_data = daily.joinpath("netcdf")
    output_data.mkdir(parents=True, exist_ok=True)
    merged = daily.joinpath("merged")
    merged.mkdir(parents=True, exist_ok=True)
    final = daily.joinpath("final")
    final.mkdir(parents=True, exist_ok=True)

    convert_flat_files(
        source_files=origin_files,
        output_folder=output_data,
        variables=var_codes,
        mode=time_step,
        n_workers=n_workers,
    )

    merge_converted_variables(
        source_files=output_data,
        output_folder=merged,
        variables=var_codes,
        n_workers=n_workers,
    )

    aggregate_stations(
        source_files=merged,
        output_folder=final,
        time_step=time_step,
        variables=var_codes,
        mf_dataset_freq="10YS",
        temp_directory=daily,
        n_workers=n_workers,
    )
