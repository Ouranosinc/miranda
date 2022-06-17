from os import getenv
from pathlib import Path

from miranda.eccc import (
    aggregate_stations,
    convert_flat_files,
    merge_converted_variables,
)

# from functools import partial
# from multiprocessing import Pool
# from miranda.eccc import convert_daily_flat_files

if __name__ == "__main__":

    time_step = "daily"
    var_codes = [
        1,
        2,
        3,
        4,
        5,
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

    station_file = source_data.joinpath("swob-xml_station_list.csv")
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
        n_workers=1,
    )

    merge_converted_variables(source=output_data, destination=merged)

    aggregate_stations(
        source_files=merged,
        output_folder=final,
        variables=var_codes,
        station_metadata=station_file,
        time_step=time_step,
        mf_dataset_freq="10YS",
        temp_directory=daily,
    )
