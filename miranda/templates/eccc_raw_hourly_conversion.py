from pathlib import Path

from miranda.eccc import aggregate_stations
from miranda.eccc import convert_hourly_flat_files
from miranda.eccc import merge_converted_variables

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
        123,
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
    # station_file = "/media/sf_VMshare/Trevor/data/Station Inventory EN.csv"
    # source_data = Path("/home/travis/doris_home/logan/scen3/smith/eccc")
    source_data = Path("/scen3/smith/eccc_converted")

    station_file = source_data.joinpath("Station Inventory EN.csv")
    origin_files = source_data.parent.joinpath("eccc_source/20200618")

    hourly = source_data.joinpath("hourly")
    output_data = hourly.joinpath("netcdf")
    output_data.mkdir(parents=True, exist_ok=True)
    merged = hourly.joinpath("merged")
    merged.mkdir(parents=True, exist_ok=True)
    final = hourly.joinpath("final")
    final.mkdir(parents=True, exist_ok=True)

    convert_hourly_flat_files(
        source_files=origin_files, output_folder=output_data, variables=var_codes
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
