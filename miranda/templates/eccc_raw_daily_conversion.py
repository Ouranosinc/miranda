from pathlib import Path

from miranda.eccc import aggregate_stations
from miranda.eccc import convert_daily_flat_files
from miranda.eccc import merge_converted_variables

if __name__ == "__main__":

    time_step = "daily"
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
    # station_file = "/media/sf_VMshare/Trevor/data/Station Inventory EN.csv"
    # source_data = Path("/home/travis/doris_home/logan/scen3/smith/eccc")
    source_data = Path("/home/tjs/Desktop/ec_data/ec")
    station_file = source_data.joinpath("Station Inventory EN.csv")
    origin_files = source_data.joinpath("source")

    daily = source_data.joinpath("daily")
    output_data = daily.joinpath("netcdf")
    output_data.mkdir(parents=True, exist_ok=True)
    merged = daily.joinpath("merged")
    merged.mkdir(parents=True, exist_ok=True)
    final = daily.joinpath("final")
    final.mkdir(parents=True, exist_ok=True)

    convert_daily_flat_files(
        source_files=source_data, output_folder=source_data, variables=var_codes
    )

    merge_converted_variables(source=output_data, destination=merged)

    aggregate_stations(
        source_files=merged,
        output_folder=final,
        variables=var_codes,
        station_metadata=station_file,
        time_step="daily",
        mf_dataset_freq="10YS",
    )
