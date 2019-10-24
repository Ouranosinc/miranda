from miranda.eccc import aggregate_hourly_nc_files
from miranda.eccc import convert_hourly_ec_files

if __name__ == "__main__":
    # indique quelle valeurs est gardee lorsqu'on a des fichiers avec plusieurs
    # valeurs valides pour les memes eccc/dates
    keep_double = "first"
    # keep_double = "last"

    var_name = "hourly_rainfall"  # "dry_bulb_temperature"  # "precipitation_amount"
    station_file = "/home/tjs/Desktop/ec_data/Station Inventory EN.csv"
    source_data = "/home/tjs/Desktop/ec_data/eccc_all"

    aggregate_hourly_nc_files(
        source_files=source_data,
        output_file=source_data,
        variable_name=var_name,
        double_handling=keep_double,
        station_inventory=station_file,
    )

    convert_hourly_ec_files(
        source_files=source_data, output_folder=source_data, variable_name=var_name
    )
