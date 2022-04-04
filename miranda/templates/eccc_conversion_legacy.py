from datetime import date
from pathlib import Path

from miranda.eccc import aggregate_stations, convert_hourly_flat_files

if __name__ == "__main__":

    var_names = [
        "atmospheric_pressure",
        "wind_speed",
        "relative_humidity",
        "dry_bulb_temperature",
        "freezing_rain",
        "ice_pellet_presence",
        "rainfall_amount",
        "precipitation_flux",
    ]
    station_file = "/home/tjs/Desktop/ec_data/Station Inventory EN.csv"
    source_data = Path("/home/tjs/Desktop/ec_data/eccc_all")

    convert_hourly_flat_files(
        source_files=source_data, output_folder=source_data, variables=var_names
    )

    for var in var_names:
        out_file = source_data.joinpath(
            "{}_eccc_hourly_{}".format(var, date.today().strftime("%Y%m%d"))
        )
        aggregate_stations(
            source_files=source_data,
            output_folder=source_data.joinpath("output"),
            variables=var,
            station_metadata=station_file,
        )
