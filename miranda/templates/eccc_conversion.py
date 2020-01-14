from datetime import date
from pathlib import Path

from miranda.eccc import aggregate_nc_files
from miranda.eccc import convert_flat_files

if __name__ == "__main__":
    # indique quelle valeurs est gardee lorsqu'on a des fichiers avec plusieurs
    # valeurs valides pour les memes eccc/dates
    # keep_double = "first"
    # keep_double = "last"

    # var_names = [
    #     "hourly_rainfall",
    #     "dry_bulb_temperature",
    #     "precipitation_amount",
    #     "freezing_rain",
    #     "ice_pellet_presence",
    # ]
    var_names = [
        "station_pressure",
        "wind_speed",
        "relative_humidity",
        "dry_bulb_temperature",
        "freezing_rain",
        "ice_pellet_presence",
        "hourly_rainfall",
        "precipitation_amount",
    ]
    station_file = "/home/tjs/Desktop/ec_data/Station Inventory EN.csv"
    source_data = Path("/home/tjs/Desktop/ec_data/eccc_all")

    convert_flat_files(
        source_files=source_data, output_folder=source_data, variables=var_names
    )

    for var in var_names:
        out_file = source_data.joinpath(
            "{}_eccc_hourly_{}".format(var, date.today().strftime("%Y%m%d"))
        )
        aggregate_nc_files(
            source_files=source_data,
            output_file=out_file,
            variables=var,
            station_inventory=station_file,
        )
