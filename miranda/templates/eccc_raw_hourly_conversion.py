from os import getenv
from pathlib import Path

from miranda.eccc import (
    aggregate_stations,
    convert_hourly_flat_files,
    merge_converted_variables,
)

# from dask.distributed import Client


if __name__ == "__main__":

    time_step = "hourly"
    var_codes = [
        76,
        77,
        78,
        79,
        80,
        89,
        94,
        107,
        108,
        109,
        110,
        123,
        133,
        156,
        262,
        263,
        264,
        265,
        266,
        267,
        268,
        269,
        270,
        271,
        272,
        273,
        274,
        275,
        276,
        277,
        278,
        279,
        280,
    ]

    var_codes = [
        262,
        263,
        264,
        265,
        266,
        267,
        268,
        269,
        270,
        271,
    ]

    in_files = getenv("in")
    source_data = Path(in_files)

    station_file = source_data.joinpath("swob-xml_station_list.csv")
    origin_files = source_data.joinpath("source")

    hourly = source_data.joinpath("hourly")
    output_data = hourly.joinpath("netcdf")
    output_data.mkdir(parents=True, exist_ok=True)
    merged = hourly.joinpath("merged")
    merged.mkdir(parents=True, exist_ok=True)
    final = hourly.joinpath("final")
    final.mkdir(parents=True, exist_ok=True)

    convert_hourly_flat_files(
        source_files=origin_files,
        output_folder=output_data,
        variables=var_codes,
        n_workers=1,
    )

    merge_converted_variables(source=output_data, destination=merged)

    aggregate_stations(
        source_files=merged,
        output_folder=final,
        variables=var_codes,
        station_metadata=station_file,
        time_step="hourly",
        mf_dataset_freq="5YS",
        temp_directory=hourly,
    )
