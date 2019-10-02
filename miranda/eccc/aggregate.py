######################################################################
# S.Biner, Ouranos, mai 2019
#
# methodologie
#
# on rassemble les fichiers netcdf des differentes eccc en un
# seul fichier netCDF.
#
# obtenu via http://climate.weather.gc.ca/index_e.html en cliquant sur 'about the data'
#######################################################################
import logging.config
import time
from datetime import datetime as dt
from pathlib import Path
from typing import Union

import cftime
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

from miranda.scripting import LOGGING_CONFIG
from miranda.utils import eccc_hourly_variable_metadata

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["aggregate_hourly_nc_files"]


def aggregate_hourly_nc_files(
    source_files: Union[str, Path],
    output_file: Union[str, Path],
    variable_name: str,
    station_inventory: Union[str, Path] = None,
    include_flags: bool = True,
    double_handling: str = "first",
):
    func_time = time.time()
    info = eccc_hourly_variable_metadata(variable_name)

    if not station_inventory:
        raise RuntimeError(
            "Download the data from ECCC's Google Drive at:\n"
            "https://drive.google.com/open?id=1egfzGgzUb0RFu_EE5AYFZtsyXPfZ11y2"
        )
    if double_handling not in ["first", "last"]:
        raise ValueError

    # On dresse la liste des eccc pour lesquelles on a des metadonnees
    df_inv = pd.read_csv(station_inventory, header=3)
    station_inventory = list(df_inv["Climate ID"].values)
    nc_data_source = Path(source_files).joinpath("split")

    # On limite le travail aux donnees et metadonnees disponibles
    rep_nc = nc_data_source.joinpath(variable_name).rglob("*.nc")
    station_file_codes = {f.name.split("_")[0] for f in rep_nc}
    stations_to_keep = set(station_file_codes).intersection(set(station_inventory))
    rejected_stations = set(station_file_codes).difference(set(station_inventory))
    valid_stations = list(sorted(stations_to_keep))
    valid_stations_count = len(valid_stations)
    logging.warning(
        "Files exist for {} ECCC stations. Metadata found for {} stations. Rejecting {} stations.".format(
            len(station_file_codes), valid_stations_count, len(rejected_stations)
        )
    )
    logging.warning(
        "Rejected station codes are the following: {}.".format(
            ", ".join(rejected_stations)
        )
    )

    # On va chercher les limites temporelles des fichiers
    list_files_to_combine = []
    for s in valid_stations:
        files = Path(nc_data_source).rglob("{}_{}*.nc".format(s, variable_name))
        files = [f.name for f in files]
        list_files_to_combine += files
    list_start_years = [int(f.split("_")[-2]) for f in list_files_to_combine]
    list_end_years = [int(f.split("_")[-1][:4]) for f in list_files_to_combine]
    year_start = np.min(list_start_years)
    year_end = np.max(list_end_years)

    # calcul des dimensions du fichier de sortie
    time_index = pd.date_range(
        start="{}-01-01".format(year_start),
        end="{}-01-01".format(year_end + 1),
        freq="H",
    )[:-1]

    # preparation du fichier de sortie
    logging.info("Preparing the NetCDF.")
    logging.info(
        "Number of eccc: {}, time steps: {}.".format(
            valid_stations_count, time_index.size
        )
    )

    # Copie/adaptation de ce qui est fait dans ec_netcdf.create_station_netcdf_file
    # des scripts utilises par Bruno Fang (et qui viennent de Blaise je pense)
    file_out = Path(output_file).joinpath(
        "{}_hourly_{}.nc".format(variable_name, double_handling)
    )
    if file_out.exists():
        file_out.unlink()
    ds = nc.Dataset(file_out, "w", format="NETCDF4")
    ds.createDimension("time", None)
    ds.createDimension("station", valid_stations_count)

    # creation de la variable
    least_significant_digit = info["least_significant_digit"]
    nc_var = ds.createVariable(
        variable_name,
        "f4",
        dimensions=("time", "station"),
        zlib=True,
        least_significant_digit=least_significant_digit,
        chunksizes=(100000, 1),
    )
    nc_var.units = info["unites"]

    if include_flags:
        # creation de la variable flag
        nc_flag = ds.createVariable(
            "flag",
            "c",
            dimensions=("time", "station"),
            zlib=True,
            chunksizes=(100000, 1),
        )
        nc_flag.long_name = "data flag"
        nc_flag.note = "See ECCC technical documentation for details"

    # Creation de la variable time
    nc_time = ds.createVariable("time", "f8", dimensions="time")
    tunits = time_index[0].strftime("days since %Y-%m-%d %H:%M:%S")
    nc_time.units = tunits
    time_num = cftime.date2num(time_index.to_pydatetime(), tunits)
    nc_time[:] = time_num
    ds.sync()

    # Creation des variables/coords des metadonnees des eccc
    nc_name = ds.createVariable("name", "str", dimensions="station")
    nc_lon = ds.createVariable("lon", "f8", dimensions="station")
    nc_lon.units = "decimal degrees"
    nc_lat = ds.createVariable("lat", "f8", dimensions="station")
    nc_lat.units = "decimal degrees"
    nc_province = ds.createVariable("province", "str", dimensions="station")
    nc_elevation = ds.createVariable("elevation", "f", dimensions="station")
    nc_elevation.units = "m"
    nc_climate_id = ds.createVariable("climate_id", "str", dimensions="station")
    nc_station_id = ds.createVariable("station_id", "str", dimensions="station")
    nc_wmo_id = ds.createVariable("wmo_id", "str", dimensions="station")
    nc_tc_id = ds.createVariable("tc_id", "str", dimensions="station")
    ds.sync()

    # On utilise les DataFrame pour aligner les donnees
    df_tot = pd.DataFrame(index=time_index)

    # On remplit les donnees pour chaques eccc
    for iter_station, station in enumerate(valid_stations):
        t00 = time.time()
        logging.info(
            "Collecting data from station {}/{}: Station ID {}.".format(
                iter_station + 1, valid_stations_count, station
            )
        )

        # On lit les donnees et les metadata pour la station
        df_stat = df_inv.loc[df_inv["Climate ID"] == station]
        # On s'assure de ne trouver qu'une metadonnee par station
        if df_stat.shape[0] != 1:
            logging.warning(
                "Duplicate metadata found for station {}. Skipping.".format(station)
            )
            continue
        df_stat = df_stat.iloc[0]

        # On ouvre les fichiers sources avec un traitement special pour les fichiers multi-annuels
        list_station_files = [
            f for f in nc_data_source.rglob("{}_{}*.nc".format(station, variable_name))
        ]
        single_year_files = []
        multi_year_files = []
        for file in list_station_files:
            a1 = file.name.split("_")[-2]
            a2 = file.name.split("_")[-1][:4]
            if a1 != a2:
                multi_year_files.append(file)
            else:
                single_year_files.append(file)

        # traitement des fichiers annuels
        if len(single_year_files) > 0:
            ds_single = xr.open_mfdataset(single_year_files, combine="by_coords")
            df_tot[variable_name] = ds_single[variable_name].to_dataframe()

            if include_flags:
                df_tot["flag"] = ds_single["flag"].to_dataframe()

        # ajout des fichiers multi annuels
        for fichier in multi_year_files:
            logging.info(
                "Special handling for multi-year data file: {}".format(fichier)
            )
            # Lecture des donnees et conversion en dataframe
            df = xr.open_dataset(fichier).to_dataframe()

            # On ne garde que les donnees valides
            df_valid = df.loc[~df[variable_name].isnull()]

            # Traitement special pour les donnees en double
            #
            # Dans les fichiers multi-annuels, il arrive qu'il y ait plusieurs
            # entrees pour les memes dates/eccc
            #
            # *** CA NE DEVRAIT PAS ARRIVER ***
            #
            if not df_valid.index.is_unique:
                logging.warning("Donnees en double ... on applique une patch ...")
                # Pour les temps ou on a plus d'une donnee valide, on garde la derniere
                df_valid = df_valid[not df_valid.index.duplicated(keep=double_handling)]
                assert df_valid.index.is_unique

            if df_valid.index.size != 0:
                df_tot.loc[df_valid.index, variable_name] = df_valid[variable_name]
                if include_flags:
                    df_tot.loc[df_valid.index, "flag"] = df_valid["flag"]

        # on copie les donnees dans le fichier de sortie
        nc_var[:, iter_station] = df_tot[variable_name].values

        if include_flags:
            nc_flag[:, iter_station] = df_tot["flag"].values

        logging.info("Added {:} data entries.".format(df_tot[variable_name].count()))

        # on enleve les donnees de df_tot
        df_tot.drop(columns=[variable_name])
        if include_flags:
            df_tot.drop(columns=["flag"])

        # on ajoute les meta donnees
        nc_name[iter_station] = df_stat["Name"]
        nc_province[iter_station] = df_stat["Province"]
        nc_lat[iter_station] = df_stat["Latitude (Decimal Degrees)"]
        nc_lon[iter_station] = df_stat["Longitude (Decimal Degrees)"]
        nc_elevation[iter_station] = df_stat["Elevation (m)"]
        nc_climate_id[iter_station] = str(df_stat["Climate ID"])
        nc_station_id[iter_station] = str(df_stat["Station ID"])
        nc_wmo_id[iter_station] = str(df_stat["WMO ID"])
        nc_tc_id[iter_station] = str(df_stat["TC ID"])
        ds.sync()

        logging.info("Transaction duration: {:.2f} seconds.".format(time.time() - t00))

    ds.Conventions = "CF-1.5"

    ds.title = "Environment and Climate Change Canada (ECCC) weather eccc"
    ds.history = "{}: Merged from multiple individual station files to n-dimensional array.".format(
        dt.now().strftime("%Y-%m-%d %X")
    )
    ds.version = "v{}".format(dt.now().strftime("%Y.%m"))
    ds.institution = "Environment and Climate Change Canada (ECCC)"
    ds.source = (
        "Weather Station data <ec.services.climatiques-climate.services.ec@canada.ca>"
    )
    ds.references = "https://climate.weather.gc.ca/doc/Technical_Documentation.pdf"
    ds.comment = "Acquired on demand from data specialists at ECCC Climate Services / Services Climatiques"
    ds.redistribution = "Redistribution policy unknown. For internal use only."

    ds.close()
    logging.warning(
        "Process completed in {:.2f} seconds".format(time.time() - func_time)
    )


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
