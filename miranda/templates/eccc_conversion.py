from miranda.eccc import aggregate_nc_files
from miranda.eccc import convert_flat_files

if __name__ == "__main__":
    # indique quelle valeurs est gardee lorsqu'on a des fichiers avec plusieurs
    # valeurs valides pour les memes eccc/dates
    keep_double = "first"
    # keep_double = "last"

    var_names = ["hourly_rainfall", "dry_bulb_temperature", "precipitation_amount"]
    station_file = "/home/tjs/Desktop/ec_data/Station Inventory EN.csv"
    source_data = "/home/tjs/Desktop/ec_data/eccc_all"

    convert_flat_files(
        source_files=source_data, output_folder=source_data, variables=var_names
    )

    aggregate_nc_files(
        source_files=source_data,
        output_file=source_data,
        variables=var_names,
        double_handling=keep_double,
        station_inventory=station_file,
    )
