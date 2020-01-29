from datetime import date
from pathlib import Path

from miranda.eccc import aggregate_nc_files
from miranda.eccc import convert_hourly_flat_files
from miranda.utils import eccc_cf_hourly_metadata

if __name__ == "__main__":

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
    station_file = "/home/tjs/Desktop/ec_data/Station Inventory EN.csv"
    source_data = Path("/home/tjs/Desktop/ec_data/eccc_all")

    convert_hourly_flat_files(
        source_files=source_data, output_folder=source_data, variables=var_codes
    )

    for var in var_codes:
        var_name = eccc_cf_hourly_metadata(var)["nc_name"]
        out_file = source_data.joinpath(
            "{}_eccc_hourly_{}".format(var_name, date.today().strftime("%Y%m%d"))
        )
        aggregate_nc_files(
            source_files=source_data,
            output_file=out_file,
            variables=var,
            station_inventory=station_file,
        )
