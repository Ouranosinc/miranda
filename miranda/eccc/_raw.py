######################################################################
# S.Biner, Ouranos, mai 2019
#
# methodologie
#
# 1) on rassemble les fichiers netcdf des differentes eccc en un seul fichier netCDF.
#
# 2) on scan les fichiers sources annuels en cherchant une variable et on sauve
# ce qu'on trouve dans des fichiers netcdf. On applique aussi les flags
# et on fait les changements d'unites
#
# obtenu via http://climate.weather.gc.ca/index_e.html en cliquant sur 'about the data'
#######################################################################
import logging
import time
from calendar import monthrange
from datetime import datetime as dt
from logging import config
from pathlib import Path
from typing import List
from typing import Union

import cftime
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

from miranda.scripting import LOGGING_CONFIG
from miranda.utils import eccc_cf_daily_metadata
from miranda.utils import eccc_hourly_variable_metadata

config.dictConfig(LOGGING_CONFIG)
__all__ = [
    "aggregate_nc_files",
    "convert_hourly_flat_files",
    "convert_daily_flat_files",
]


def convert_hourly_flat_files(
    source_files: Union[str, Path],
    output_folder: Union[str, Path],
    variables: Union[str, List[str]],
    missing_value: int = -9999,
) -> None:
    """

    Parameters
    ----------
    source_files : Union[str, Path]
    output_folder : Union[str, Path]
    variables : Union[str, List[str]
    missing_value : int

    Returns
    -------
    None
    """
    func_time = time.time()

    if isinstance(variables, str):
        variables = [variables]

    for variable_name in variables:
        info = eccc_hourly_variable_metadata(variable_name)
        variable_code = info["code_var"]

        # preparing the data extraction
        titre_colonnes = "code year month day code_var ".split()
        for i in range(1, 25):
            titre_colonnes.append("D{:0n}".format(i))
            titre_colonnes.append("F{:0n}".format(i))

        rep_nc = Path(output_folder).joinpath(variable_name)
        rep_nc.mkdir(parents=True, exist_ok=True)

        # Loop on the files
        if 262 < int(variable_code) <= 280:
            list_files = Path(source_files).rglob("HLY*RCS*.gz")
        else:
            list_files = Path(source_files).rglob("HLY*.gz")

        for fichier in list_files:
            logging.info("Processing file: {}.".format(fichier))

            # Create a dataframe from the files
            df = pd.read_fwf(
                fichier, widths=[7, 4, 2, 2, 3] + [6, 1] * 24, names=titre_colonnes
            )

            # Loop through the station codes
            l_codes = df["code"].unique()
            for code in l_codes:
                df_code = df[df["code"] == code]

                # on arrete si la variable n'est pas presente
                if variable_code not in df_code["code_var"].unique():
                    continue

                # on fait le traitement
                logging.info(
                    "Converting {} for station code: {}".format(variable_name, code)
                )

                # on va chercher le dataframe de la variable
                df_var = df_code[df_code["code_var"] == variable_code].copy()

                # application du masque selon la valeur manquante
                df_var = df_var.replace(missing_value, np.nan)

                # decodage des valeurs et flags
                dfd = df_var.loc[:, ["D{:0n}".format(i) for i in range(1, 25)]]
                dff = df_var.loc[:, ["F{:0n}".format(i) for i in range(1, 25)]]

                # enleve flag == nan
                dff = dff.fillna("")

                # utilise flag pour masquer valeurs
                val = np.asfarray(dfd.values)
                flag = dff.values
                mask = np.isin(flag, info["missing_flags"])
                val[mask] = np.nan

                # traitement des unites
                val = val * info["scale_factor"] + info["add_offset"]

                # Create the DataArray
                dates = dict(time=list())
                for index, row in df_var.iterrows():
                    for h in range(0, 24):
                        dates["list"].append(
                            dt(int(row.year), int(row.month), int(row.day), h)
                        )

                ds = xr.Dataset()
                da_val = xr.DataArray(val.flatten(), coords=dates, dims=["time"])
                da_val = da_val.rename(variable_name)
                da_val.attrs["unites"] = info["unites"]
                da_val.attrs["id"] = code
                da_val.attrs["element_number"] = variable_code

                da_flag = xr.DataArray(flag.flatten(), coords=dates, dims=["time"])
                ds[variable_name] = da_val
                ds["flag"] = da_flag

                # save the file in NetCDF format
                start_year = ds.time.dt.year.values[0]
                end_year = ds.time.dt.year.values[-1]
                if start_year == end_year:
                    f_nc = "{c}_{vc}_{v}_{sy}.nc".format(
                        c=code, vc=variable_code, v=variable_name, sy=start_year
                    )
                else:
                    f_nc = "{c}_{vc}_{v}_{sy}_{ey}.nc".format(
                        c=code,
                        vc=variable_code,
                        v=variable_name,
                        sy=start_year,
                        ey=end_year,
                    )

                ds.attrs["Conventions"] = "CF-1.5"

                ds.attrs[
                    "title"
                ] = "Environment and Climate Change Canada (ECCC) weather eccc"
                ds.attrs[
                    "history"
                ] = "{}: Merged from multiple individual station files to n-dimensional array.".format(
                    dt.now().strftime("%Y-%m-%d %X")
                )
                ds.attrs["version"] = "v{}".format(dt.now().strftime("%Y.%m"))
                ds.attrs["institution"] = "Environment and Climate Change Canada (ECCC)"
                ds.attrs[
                    "source"
                ] = "Weather Station data <ec.services.climatiques-climate.services.ec@canada.ca>"
                ds.attrs[
                    "references"
                ] = "https://climate.weather.gc.ca/doc/Technical_Documentation.pdf"
                ds.attrs[
                    "comment"
                ] = "Acquired on demand from data specialists at ECCC Climate Services / Services Climatiques"
                ds.attrs[
                    "redistribution"
                ] = "Redistribution policy unknown. For internal use only."

                ds.to_netcdf(rep_nc.joinpath(f_nc))

    logging.warning(
        "Process completed in {:.2f} seconds".format(time.time() - func_time)
    )


def convert_daily_flat_files(
    source_files: Union[str, Path],
    output_folder: Union[str, Path],
    variables: Union[str, int, List[Union[str, int]]],
    missing_value: int = -9999,
) -> None:
    """

    Parameters
    ----------
    source_files : Union[str, Path]
    output_folder : Union[str, Path]
    variables : Union[str, List[str]
      Variable codes (001, 002, 103, etc.)
    missing_value : int

    Returns
    -------
    None
    """
    func_time = time.time()

    if isinstance(variables, str):
        variables = [variables]

    for variable_code in variables:
        info = eccc_cf_daily_metadata(variable_code)
        variable_code = str(variable_code).zfill(3)
        variable_name = info["standard_name"]

        # on prepare l'extraction des donnees
        titre_colonnes = "code year month code_var".split()
        for i in range(1, 32):
            titre_colonnes.append("D{:0n}".format(i))
            titre_colonnes.append("F{:0n}".format(i))

        # preparation du repertoire de sortie
        rep_nc = Path(output_folder).joinpath(variable_name)
        rep_nc.mkdir(parents=True, exist_ok=True)

        # boucle sur les fichiers
        list_files = list()

        if isinstance(source_files, list) or Path(source_files).is_file():
            list_files.append(source_files)
        else:
            list_files.extend(Path(source_files).rglob("*DLY*"))

        for fichier in list_files:
            logging.info("Processing file: {}.".format(fichier))

            # Create a dataframe from the files
            df = pd.read_fwf(
                fichier,
                widths=[7, 4, 2, 3] + [6, 1] * 31,
                names=titre_colonnes,
                dtype={"year": int, "month": int, "code_var": str},
            )

            # Loop through the station codes
            l_codes = df["code"].unique()
            for code in l_codes:
                df_code = df[df["code"] == code]

                # on arrete si la variable n'est pas presente
                if variable_code not in df_code["code_var"].unique():
                    continue

                # on fait le traitement
                logging.info(
                    "Converting {} for station code: {}".format(variable_name, code)
                )

                # on va chercher le dataframe de la variable
                df_var = df_code[df_code["code_var"] == variable_code].copy()

                # application du masque selon la valeur manquante
                df_var = df_var.replace(missing_value, np.nan)

                # Decoding the values and flags
                dfd = df_var.loc[:, ["D{:0n}".format(i) for i in range(1, 32)]]
                dff = df_var.loc[:, ["F{:0n}".format(i) for i in range(1, 32)]]

                # enleve flag == nan
                dff = dff.fillna("")

                # utilise flag pour masquer valeurs
                val = np.asfarray(dfd.values)
                flag = dff.values
                mask = np.isin(flag, info["missing_flags"])
                val[mask] = np.nan

                # Adjust units
                val = val * info["scale_factor"] + info["add_offset"]

                # Create the DataArray and concatenate values and flags based on day-length of months
                date_range = dict(time=list())
                value_days = list()
                flag_days = list()
                for i, (index, row) in enumerate(df_var.iterrows()):
                    period = pd.Period(year=row.year, month=row.month, freq="M")
                    dates = pd.Series(
                        pd.date_range(
                            start=period.start_time, end=period.end_time, freq="D"
                        )
                    )
                    date_range["time"].extend(dates)

                    value_days.extend(val[i][range(monthrange(row.year, row.month)[1])])
                    flag_days.extend(flag[i][range(monthrange(row.year, row.month)[1])])

                ds = xr.Dataset()
                da_val = xr.DataArray(value_days, coords=date_range, dims=["time"])
                da_val = da_val.rename(variable_name)
                da_val.attrs["units"] = info["nc_units"]
                da_val.attrs["id"] = code
                da_val.attrs["element_number"] = variable_code

                da_flag = xr.DataArray(flag_days, coords=date_range, dims=["time"])
                ds[variable_name] = da_val
                ds["flag"] = da_flag

                # save as a NetCDF file
                start_year = ds.time.dt.year.values[0]
                end_year = ds.time.dt.year.values[-1]
                if start_year == end_year:
                    f_nc = "{c}_{vc}_{v}_{sy}.nc".format(
                        c=code, vc=variable_code, v=variable_name, sy=start_year
                    )
                else:
                    f_nc = "{c}_{vc}_{v}_{sy}_{ey}.nc".format(
                        c=code,
                        vc=variable_code,
                        v=variable_name,
                        sy=start_year,
                        ey=end_year,
                    )

                ds.attrs["Conventions"] = "CF-1.7"

                ds.attrs[
                    "title"
                ] = "Environment and Climate Change Canada (ECCC) weather eccc"
                ds.attrs[
                    "history"
                ] = "{}: Merged from multiple individual station files to n-dimensional array.".format(
                    dt.now().strftime("%Y-%m-%d %X")
                )
                ds.attrs["version"] = "v{}".format(dt.now().strftime("%Y.%m"))
                ds.attrs["institution"] = "Environment and Climate Change Canada (ECCC)"
                ds.attrs[
                    "source"
                ] = "Weather Station data <ec.services.climatiques-climate.services.ec@canada.ca>"
                ds.attrs[
                    "references"
                ] = "https://climate.weather.gc.ca/doc/Technical_Documentation.pdf"
                ds.attrs[
                    "comment"
                ] = "Acquired on demand from data specialists at ECCC Climate Services / Services Climatiques"
                ds.attrs[
                    "redistribution"
                ] = "Redistribution policy unknown. For internal use only."

                ds.to_netcdf(rep_nc.joinpath(f_nc))

    logging.warning(
        "Process completed in {:.2f} seconds".format(time.time() - func_time)
    )


# TODO: Adjust this function to allow for hourly and daily data aggregation
def aggregate_nc_files(
    source_files: Union[str, Path],
    output_file: Union[str, Path],
    variables: Union[str, List[str]],
    time_step: str = "h",
    station_inventory: Union[str, Path] = None,
    include_flags: bool = True,
    double_handling: str = "first",
) -> None:
    """

    Parameters
    ----------
    source_files: Union[str, Path]
    output_file: Union[str, Path]
    variables: Union[str, List[str]]
    time_step: str
    station_inventory: Union[str, Path]
    include_flags: bool
    double_handling: str

    Returns
    -------
    None
    """
    func_time = time.time()

    if not station_inventory:
        raise RuntimeError(
            "Download the data from ECCC's Google Drive at:\n"
            "https://drive.google.com/open?id=1egfzGgzUb0RFu_EE5AYFZtsyXPfZ11y2"
        )
    if double_handling not in ["first", "last"]:
        raise ValueError

    if isinstance(variables, str):
        variables = [variables]

    if isinstance(source_files, str):
        source_files = Path(source_files)

    if time_step.lower() in ["h", "hour", "hourly"]:
        hourly = True
    elif time_step.lower() in ["d", "day", "daily"]:
        hourly = False
    else:
        raise ValueError("Time step must be `hourly` or `daily`.")

    for variable_name in variables:
        if hourly:
            info = eccc_hourly_variable_metadata(variable_name)
        else:
            info = eccc_cf_daily_metadata(variable_name)

        # On dresse la liste des eccc pour lesquelles on a des metadonnees
        df_inv = pd.read_csv(station_inventory, header=3)
        station_inventory = list(df_inv["Climate ID"].values)

        # On limite le travail aux donnees et metadonnees disponibles
        rep_nc = source_files.joinpath(variable_name).rglob("*.nc")
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
            files = [
                f.name
                for f in Path(source_files).rglob("{}_{}*.nc".format(s, variable_name))
            ]
            list_files_to_combine += files
        list_start_years = [int(f.split("_")[-2]) for f in list_files_to_combine]
        list_end_years = [int(f.split("_")[-1][:4]) for f in list_files_to_combine]
        year_start, year_end = np.min(list_start_years), np.max(list_end_years)

        # Calculate the dimensions of the output NetCDF
        time_index = pd.date_range(
            start="{}-01-01".format(year_start),
            end="{}-01-01".format(year_end + 1),
            freq=["H" if hourly else "D"],
        )[:-1]

        logging.info("Preparing the NetCDF.")
        logging.info(
            "Number of eccc: {}, time steps: {}.".format(
                valid_stations_count, time_index.size
            )
        )

        # Copie/adaptation de ce qui est fait dans ec_netcdf.create_station_netcdf_file
        # des scripts utilises par Bruno Fang (et qui viennent de Blaise je pense)
        if not Path(output_file).exists():
            Path(output_file).mkdir()

        file_out = Path(output_file).joinpath(
            "{}_{}_{}.nc".format(
                variable_name, "hourly" if hourly else "daily", double_handling
            )
        )
        if file_out.exists():
            file_out.unlink()

        ds = netCDF4.Dataset(file_out, "w", format="NETCDF4")
        ds.createDimension("time", None)
        ds.createDimension("station", valid_stations_count)

        # creation de la variable
        least_significant_digit = info["least_significant_digit"]
        nc_var = ds.createVariable(
            variable_name,
            "f4",  # TODO: create a mapping of the data type for floats, int and bools to nc_datatypes
            dimensions=("time", "station"),
            zlib=True,
            least_significant_digit=least_significant_digit,
            chunksizes=(100000, 1),
        )
        nc_var.units = info["nc_units"]

        nc_flag = None
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

        # Create the time variable
        nc_time = ds.createVariable("time", "f8", dimensions="time")
        tunits = time_index[0].strftime("days since %Y-%m-%d %H:%M:%S")
        nc_time.units = tunits
        time_num = cftime.date2num(time_index.to_pydatetime(), tunits)
        nc_time[:] = time_num
        ds.sync()

        # Write out the appropriate ECCC metadata
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

        # Use a Pandas DataFrame to align the data
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
                f
                for f in source_files.rglob("{}_{}*.nc".format(station, variable_name))
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
                try:
                    ds_single = xr.open_mfdataset(
                        single_year_files, combine="by_coords"
                    )
                except ValueError:
                    logging.exception(
                        "Unable to read by coords. Skipping station: {}".format(station)
                    )
                    continue
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
                    df_valid = df_valid[
                        not df_valid.index.duplicated(keep=double_handling)
                    ]
                    assert df_valid.index.is_unique

                if df_valid.index.size != 0:
                    df_tot.loc[df_valid.index, variable_name] = df_valid[variable_name]
                    if include_flags:
                        df_tot.loc[df_valid.index, "flag"] = df_valid["flag"]

            # on copie les donnees dans le fichier de sortie
            nc_var[:, iter_station] = df_tot[variable_name].values

            if include_flags:
                nc_flag[:, iter_station] = df_tot["flag"].values

            logging.info(
                "Added {:} data entries.".format(df_tot[variable_name].count())
            )

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

            logging.info(
                "Transaction duration: {:.2f} seconds.".format(time.time() - t00)
            )

        ds.Conventions = "CF-1.7"

        ds.title = "Environment and Climate Change Canada (ECCC) weather eccc"
        ds.history = "{}: Merged from multiple individual station files to n-dimensional array.".format(
            dt.now().strftime("%Y-%m-%d %X")
        )
        ds.version = "v{}".format(dt.now().strftime("%Y.%m"))
        ds.institution = "Environment and Climate Change Canada (ECCC)"
        ds.source = "Weather Station data <ec.services.climatiques-climate.services.ec@canada.ca>"
        ds.references = "https://climate.weather.gc.ca/doc/Technical_Documentation.pdf"
        ds.comment = "Acquired on demand from data specialists at ECCC Climate Services / Services Climatiques"
        ds.redistribution = "Redistribution policy unknown. For internal use only."

        ds.close()
    logging.warning(
        "Process completed in {:.2f} seconds".format(time.time() - func_time)
    )
