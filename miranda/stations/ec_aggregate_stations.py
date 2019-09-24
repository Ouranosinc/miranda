import logging.config
import time
from pathlib import Path
from typing import Union

import cftime
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

from miranda.utils import get_info_var

# from miranda.scripting import LOGGING_CONFIG
# logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["aggregate_hourly_nc_files"]


def aggregate_hourly_nc_files(
    source_files: Union[str, Path],
    output_file: Union[str, Path],
    variable_name: str,
    station_inventory: Union[str, Path] = None,
    include_flags: bool = True,
    double_handling: str = "first",
):
    info = get_info_var(variable_name)

    if not station_inventory:
        raise RuntimeError(
            "Download the data from ECCC's Google Drive at:\n"
            "https://drive.google.com/open?id=1egfzGgzUb0RFu_EE5AYFZtsyXPfZ11y2"
        )
    if double_handling not in ["first", "last"]:
        raise ValueError

    # on dresse la liste des stations pour lesquelles on a des metadonnees
    df_inv = pd.read_csv(station_inventory, header=3)
    station_inventory = list(df_inv["Climate ID"].values)
    nc_data_source = Path(source_files).joinpath("split")

    # on limite le travail aux donnees et metadonnees disponibles
    rep_nc = nc_data_source.joinpath(nom_var).rglob("*.nc")
    list_of_stations = {f.name.split("_")[0] for f in rep_nc}
    stations_to_keep = set(list_of_stations).intersection(set(station_inventory))
    rejected_stations = set(list_of_stations).difference(set(station_inventory))
    list_of_stations = list(sorted(stations_to_keep))
    logging.warning(
        "on garde {} stations et on rejette {} stations".format(
            len(list_of_stations), len(rejected_stations)
        )
    )

    # on va chercher les limites temporelles des fichiers
    list_files_to_combine = []
    for s in list_of_stations:
        files = Path(nc_data_source).rglob("{}_{}*.nc".format(s, nom_var))
        files = [f.name for f in files]
        list_files_to_combine += files
    list_start_years = [int(f.split("_")[-2]) for f in list_files_to_combine]
    list_end_years = [int(f.split("_")[-1][:4]) for f in list_files_to_combine]
    year_start = np.min(list_start_years)
    year_end = np.max(list_end_years)

    # calcul des dimensions du fichier de sortie
    nb_stations = len(list_of_stations)
    time_index = pd.date_range(
        start="{}-01-01".format(year_start),
        end="{}-01-01".format(year_end + 1),
        freq="H",
    )[:-1]

    # preparation du fichier de sortie
    logging.info("Preparation du fichier de sortie.")
    logging.info("Dimensions station:{}, time:{}".format(nb_stations, time_index.size))

    # copie/adaptation de ce qui est fait dans ec_netcdf.create_station_netcdf_file
    # des scripts utilises par Bruno Fang (et qui viennent de Blaise je pense)
    file_out = Path(output_file).joinpath(
        "{}_hourly_{}.nc".format(nom_var, double_handling)
    )
    if file_out.exists():
        file_out.unlink()
    dso = nc.Dataset(file_out, "w", format="NETCDF4")
    dso.createDimension("time", None)
    dso.createDimension("station", nb_stations)

    # creation de la variable
    least_significant_digit = info["least_significant_digit"]
    nc_var = dso.createVariable(
        nom_var,
        "f4",
        dimensions=("time", "station"),
        zlib=True,
        least_significant_digit=least_significant_digit,
        chunksizes=(100000, 1),
    )
    nc_var.units = info["unites"]

    if include_flags:
        # creation de la variable flag
        nc_flag = dso.createVariable(
            "flag",
            "c",
            dimensions=("time", "station"),
            zlib=True,
            chunksizes=(100000, 1),
        )
        nc_flag.long_name = "data flag"
        nc_flag.note = "see ECCC technical documentation for details"

    # creation de la variable time
    nc_time = dso.createVariable("time", "f8", dimensions="time")
    tunits = time_index[0].strftime("days since %Y-%m-%d %H:%M:%S")
    nc_time.units = tunits
    time_num = cftime.date2num(time_index.to_pydatetime(), tunits)
    nc_time[:] = time_num
    dso.sync()

    # creation des variables/coords des metadonnees des stations
    nc_nom = dso.createVariable("nom", "str", dimensions="station")
    nc_lon = dso.createVariable("lon", "f8", dimensions="station")
    nc_lon.units = "decimal degrees"
    nc_lat = dso.createVariable("lat", "f8", dimensions="station")
    nc_lat.units = "decimal degrees"
    nc_province = dso.createVariable("province", "str", dimensions="station")
    nc_elevation = dso.createVariable("elevation", "f", dimensions="station")
    nc_elevation.units = "m"
    nc_climate_id = dso.createVariable("climate_id", "str", dimensions="station")
    nc_station_id = dso.createVariable("station_id", "str", dimensions="station")
    nc_wmo_id = dso.createVariable("wmo_id", "str", dimensions="station")
    nc_tc_id = dso.createVariable("tc_id", "str", dimensions="station")
    dso.sync()

    # on remplit les donnees pour chaques stations
    for n_s, station in enumerate(list_of_stations):
        t00 = time.time()
        logging.info(
            "traitement de station {}/{} : {}".format(n_s + 1, nb_stations, station)
        )

        # on lit les donnees et les metadata pour la station
        df_stat = df_inv.loc[df_inv["Climate ID"] == station]
        # on s'assure de ne trouver qu'une metadonnee par station
        assert df_stat.shape[0] == 1
        df_stat = df_stat.iloc[0]

        # on utilise les DataFrame pour aligner les donnees
        if n_s == 0:
            df_tot = pd.DataFrame(index=time_index)

        # on ouvre les fichiers sources avec un traitement special pour les
        # fichiers multi-annuels

        list_station_files = [
            f for f in nc_data_source.rglob("{}_*.nc".format(station))
        ]
        single_year_files = []
        multi_year_files = []
        for file in list_station_files:
            a1 = file.name.split("_")[2]
            a2 = file.name.split("_")[3][:4]
            if a1 != a2:
                multi_year_files.append(file)
            else:
                single_year_files.append(file)

        # traitement des fichiers annuels
        if len(single_year_files) > 0:
            ds = xr.open_mfdataset(single_year_files, combine="by_coords")
            df_tot[nom_var] = ds[nom_var].to_dataframe()

            if include_flags:
                df_tot["flag"] = ds["flag"].to_dataframe()

        # ajout des fichiers multi annuels
        for fichier in multi_year_files:
            logging.info(
                "traitement special pour fichier multi_annuels {}".format(fichier)
            )
            # lecture des donnees et conversion en dataframe
            ds = xr.open_dataset(fichier)
            df = ds.to_dataframe()

            # on ne garde que les donnees valides
            df_valid = df.loc[~df[nom_var].isnull()]

            # traitement special pour les donnees en double
            #
            # Dans les fichiers multi-annuels, il arrive qu'il y ait plusieurs
            # entrees pour les memes dates/stations
            #
            # *** CA NE DEVRAIT PAS ARRIVER ***
            #
            if not df_valid.index.is_unique:
                logging.warning("donnees en double ... on applique une patch ...")
                # pour les temps ou on a plus d'une donnee valide, on garde
                # la derniere
                df_valid = df_valid[not df_valid.index.duplicated(keep=double_handling)]
                assert df_valid.index.is_unique

            if df_valid.index.size != 0:
                df_tot.loc[df_valid.index, nom_var] = df_valid[nom_var]
                if include_flags:
                    df_tot.loc[df_valid.index, "flag"] = df_valid["flag"]

        # on copie les donnees dans le fichier de sortie
        nc_var[:, n_s] = df_tot[nom_var].values

        if include_flags:
            nc_flag[:, n_s] = df_tot["flag"].values

        logging.info("ajout de {:} donnees".format(df_tot[nom_var].count()))

        # on enleve les donnees de df_tot
        df_tot.drop(columns=[nom_var])
        if include_flags:
            df_tot.drop(columns=["flag"])

        # on ajoute les meta donnees
        nc_nom[n_s] = df_stat["Name"]
        nc_province[n_s] = df_stat["Province"]
        nc_lat[n_s] = df_stat["Latitude (Decimal Degrees)"]
        nc_lon[n_s] = df_stat["Longitude (Decimal Degrees)"]
        nc_elevation[n_s] = df_stat["Elevation (m)"]
        nc_climate_id[n_s] = str(df_stat["Climate ID"])
        nc_station_id[n_s] = str(df_stat["Station ID"])
        nc_wmo_id[n_s] = str(df_stat["WMO ID"])
        nc_tc_id[n_s] = str(df_stat["TC ID"])
        dso.sync()

        logging.info("timing: {:.2f}(s)".format(time.time() - t00))

    dso.close()


if __name__ == "__main__":
    # indique quelle valeurs est gardee lorsqu'on a des fichiers avec plusieurs
    # valeurs valides pour les memes stations/dates
    keep_double = "first"
    # keep_double = "last"

    nom_var = "tas"  # "rainfall"
    station_file = "/home/tjs/Desktop/ec_data/Station Inventory EN.csv"
    source_data = "/home/tjs/Desktop/ec_data"

    aggregate_hourly_nc_files(
        source_files=source_data,
        output_file=source_data,
        variable_name=nom_var,
        double_handling=keep_double,
        station_inventory=station_file,
    )
