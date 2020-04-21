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
from typing import Optional
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar

from miranda.scripting import LOGGING_CONFIG
from miranda.utils import eccc_cf_daily_metadata
from miranda.utils import eccc_cf_hourly_metadata

config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "aggregate_nc_files",
    "convert_hourly_flat_files",
    "convert_daily_flat_files",
]


def convert_hourly_flat_files(
    source_files: Union[str, Path],
    output_folder: Union[str, Path, List[Union[str, int]]],
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

    if isinstance(variables, (str, int)):
        variables = [variables]

    for variable_code in variables:
        info = eccc_cf_hourly_metadata(variable_code)
        variable_code = str(variable_code).zfill(3)
        variable_name = info["standard_name"]
        variable_file_name = info["nc_name"]

        # Preparing the data extraction
        titre_colonnes = "code year month day code_var ".split()
        for i in range(1, 25):
            titre_colonnes.append("D{:0n}".format(i))
            titre_colonnes.append("F{:0n}".format(i))

        rep_nc = Path(output_folder).joinpath(variable_file_name)
        rep_nc.mkdir(parents=True, exist_ok=True)

        # Loop on the files
        if 262 < int(variable_code) <= 280:
            list_files = Path(source_files).rglob("HLY*RCS*.gz")
        else:
            list_files = Path(source_files).rglob("HLY*.gz")

        errored_files = list()
        for fichier in list_files:
            logging.info("Processing file: {}.".format(fichier))

            # Create a dataframe from the files
            try:
                df = pd.read_fwf(
                    fichier,
                    widths=[7, 4, 2, 2, 3] + [6, 1] * 24,
                    names=titre_colonnes,
                    dtype={"year": int, "month": int, "day": int, "code_var": str},
                )
            except ValueError:
                logging.error(
                    "File {} was unable to be read. This is probably an issue with the file.".format(
                        fichier
                    )
                )
                errored_files.append(fichier)
                continue

            # Loop through the station codes
            l_codes = df["code"].unique()
            for code in l_codes:
                df_code = df[df["code"] == code]

                # Abort if the variable is not found
                if variable_code not in df_code["code_var"].unique():
                    logging.info(
                        "Variable `{}` not found for station code: {}. Continuing...".format(
                            variable_file_name, code
                        )
                    )
                    continue

                # Treat the data
                logging.info(
                    "Converting `{}` for station code: {}".format(
                        variable_file_name, code
                    )
                )

                # Dump the data into a DataFrame
                df_var = df_code[df_code["code_var"] == variable_code].copy()

                # Mask the data according to the missing values flag
                df_var = df_var.replace(missing_value, np.nan)

                # Decode the values and flags
                dfd = df_var.loc[:, ["D{:0n}".format(i) for i in range(1, 25)]]
                dff = df_var.loc[:, ["F{:0n}".format(i) for i in range(1, 25)]]

                # Remove the "NaN" flag
                dff = dff.fillna("")

                # Use the flag to mask the values
                val = np.asfarray(dfd.values)
                flag = dff.values
                mask = np.isin(flag, info["missing_flags"])
                val[mask] = np.nan

                # Treat according to units conversions
                val = val * info["scale_factor"] + info["add_offset"]

                # Create the DataArray
                dates = dict(time=list())
                for index, row in df_var.iterrows():
                    for h in range(0, 24):
                        dates["time"].append(
                            dt(int(row.year), int(row.month), int(row.day), h)
                        )

                ds = xr.Dataset()
                da_val = xr.DataArray(val.flatten(), coords=dates, dims=["time"])
                da_val = da_val.rename(variable_file_name)
                da_val.attrs["units"] = info["nc_units"]
                da_val.attrs["id"] = code
                da_val.attrs["element_number"] = variable_code
                da_val.attrs["standard_name"] = variable_name
                da_val.attrs["long_name"] = info["long_name"]

                da_flag = xr.DataArray(flag.flatten(), coords=dates, dims=["time"])
                da_flag.attrs["long_name"] = "data flag"
                da_flag.attrs["note"] = "See ECCC technical documentation for details"

                ds[variable_file_name] = da_val
                ds["flag"] = da_flag

                # save the file in NetCDF format
                start_year = ds.time.dt.year.values[0]
                end_year = ds.time.dt.year.values[-1]

                station_folder = rep_nc.joinpath(str(code))
                station_folder.mkdir(parents=True, exist_ok=True)

                if start_year == end_year:
                    f_nc = "{c}_{vc}_{v}_{sy}.nc".format(
                        c=code, vc=variable_code, v=variable_file_name, sy=start_year
                    )
                else:
                    f_nc = "{c}_{vc}_{v}_{sy}_{ey}.nc".format(
                        c=code,
                        vc=variable_code,
                        v=variable_file_name,
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

                ds.to_netcdf(station_folder.joinpath(f_nc))

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
    variables : Union[str, int, List[Union[str, int]]
      Variable codes (001, 002, 103, etc.)
    missing_value : int

    Returns
    -------
    None
    """
    func_time = time.time()

    if isinstance(variables, (str, int)):
        variables = [variables]

    for variable_code in variables:
        info = eccc_cf_daily_metadata(variable_code)
        variable_code = str(variable_code).zfill(3)
        nc_name = info["nc_name"]

        # Prepare the data extraction
        titre_colonnes = "code year month code_var".split()
        for i in range(1, 32):
            titre_colonnes.append("D{:0n}".format(i))
            titre_colonnes.append("F{:0n}".format(i))

        # Create the output directory
        rep_nc = Path(output_folder).joinpath(nc_name)
        rep_nc.mkdir(parents=True, exist_ok=True)

        # Loop on the files
        list_files = list()
        if isinstance(source_files, list) or Path(source_files).is_file():
            list_files.append(source_files)
        else:
            list_files.extend(Path(source_files).rglob("*DLY*"))

        errored_files = list()
        for fichier in list_files:
            logging.info("Processing file: {}.".format(fichier))

            # Create a Pandas DataFrame from the files
            try:
                df = pd.read_fwf(
                    fichier,
                    widths=[7, 4, 2, 3] + [6, 1] * 31,
                    names=titre_colonnes,
                    dtype={"year": int, "month": int, "code_var": str},
                )
            except ValueError:
                logging.error(
                    "File {} was unable to be read. This is probably an issue with the file.".format(
                        fichier
                    )
                )
                errored_files.append(fichier)
                continue

            # Loop through the station codes
            l_codes = df["code"].unique()
            for code in l_codes:
                df_code = df[df["code"] == code]

                # Abort if the variable is not present
                if variable_code not in df_code["code_var"].unique():
                    logging.info(
                        "Variable `{}` not found for station `{}` in file {}. Continuing...".format(
                            nc_name, code, fichier
                        )
                    )
                    continue

                # Perform the data treatment
                logging.info("Converting {} for station code: {}".format(nc_name, code))

                # Dump the values into a DataFrame
                df_var = df_code[df_code["code_var"] == variable_code].copy()

                # Apply the mask according to the NaN value
                df_var = df_var.replace(missing_value, np.nan)

                # Decoding the values and flags
                dfd = df_var.loc[:, ["D{:0n}".format(i) for i in range(1, 32)]]
                dff = df_var.loc[:, ["F{:0n}".format(i) for i in range(1, 32)]]

                # Remove the "NaN" flag
                dff = dff.fillna("")

                # Use the flag to mask the values
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
                da_val = da_val.rename(nc_name)
                da_val.attrs["units"] = info["nc_units"]
                da_val.attrs["id"] = code
                da_val.attrs["element_number"] = variable_code
                da_val.attrs["standard_name"] = info["standard_name"]
                da_val.attrs["long_name"] = info["long_name"]

                da_flag = xr.DataArray(flag_days, coords=date_range, dims=["time"])
                da_flag.attrs["long_name"] = "data flag"
                da_flag.attrs["note"] = "See ECCC technical documentation for details"

                ds[nc_name] = da_val
                ds["flag"] = da_flag

                # Save as a NetCDF file
                start_year = ds.time.dt.year.values[0]
                end_year = ds.time.dt.year.values[-1]

                station_folder = rep_nc.joinpath(str(code))
                station_folder.mkdir(parents=True, exist_ok=True)

                if start_year == end_year:
                    f_nc = "{c}_{vc}_{v}_{sy}.nc".format(
                        c=code, vc=variable_code, v=nc_name, sy=start_year
                    )
                else:
                    f_nc = "{c}_{vc}_{v}_{sy}_{ey}.nc".format(
                        c=code, vc=variable_code, v=nc_name, sy=start_year, ey=end_year,
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

                ds.to_netcdf(station_folder.joinpath(f_nc))

    logging.warning(
        "Process completed in {:.2f} seconds".format(time.time() - func_time)
    )


def aggregate_nc_files(
    source_files: Optional[Union[str, Path]] = None,
    output_folder: Optional[Union[str, Path]] = None,
    station_inventory: Union[str, Path] = None,
    time_step: str = "h",
    variables: Optional[Union[str, int, List[Union[str, int]]]] = None,
    include_flags: bool = True,
    mf_dataset_freq: Optional[str] = None,
) -> None:
    """

    Parameters
    ----------
    source_files: Union[str, Path]
    output_folder: Union[str, Path]
    variables: Optional[Union[str, int, List[Union[str, int]]]]
    time_step: str
    station_inventory: Union[str, Path]
    include_flags: bool
    mf_dataset_freq: Optional[str]
      Resampling frequency for creating output multi-file Datasets. E.g. 'YS': 1 year per file, '5YS': 5 y per file

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

    if isinstance(source_files, str):
        source_files = Path(source_files)

    if time_step.lower() in ["h", "hour", "hourly"]:
        hourly = True
    elif time_step.lower() in ["d", "day", "daily"]:
        hourly = False
    else:
        raise ValueError("Time step must be `h` / `hourly` or `d` / `daily`.")

    if isinstance(variables, (str, int)):
        variables = [variables]
    elif variables is None:
        if hourly:
            variables = [
                89,
                94,
                123,
            ]
            variables.extend(range(76, 81))
            variables.extend(range(262, 281))
        else:
            variables = [1, 2, 3]
            variables.extend(range(10, 26))

    for variable_code in variables:
        if hourly:
            info = eccc_cf_hourly_metadata(variable_code)
        else:
            info = eccc_cf_daily_metadata(variable_code)
        variable_name = info["nc_name"]
        logging.info(f"Merging `{variable_name}` using `{time_step}` time step.")

        # Find the ECCC stations where we have available metadata
        df_inv = pd.read_csv(station_inventory, header=3)
        station_inventory = list(df_inv["Climate ID"].values)

        # Only perform aggregation on available data with corresponding metadata
        logging.info("Performing glob and sort.")
        nclist = sorted(list(source_files.joinpath(variable_name).rglob("*.nc")))
        station_file_codes = [f.name.split("_")[0] for f in nclist]
        stations_to_keep = list(
            set(station_file_codes).intersection(set(station_inventory))
        )
        logging.info(f"{len(stations_to_keep)} stations found and sorted.")
        if len(nclist) > 0:
            ds = xr.open_mfdataset(nclist, combine="nested", concat_dim="station")
            ds = ds.assign_coords(
                station_id=xr.DataArray(station_file_codes, dims="station")
            )

            rejected_stations = set(station_file_codes).difference(
                set(station_inventory)
            )
            logging.info(f"{len(rejected_stations)} rejected due to missing metadata.")
            for r in rejected_stations:
                ds = ds.isel(station=(ds.station_id != r))
            if not include_flags:
                drop_vars = [vv for vv in ds.data_vars if "flag" in vv]
                ds = ds.drop_vars(drop_vars)

            # make sure data is in order to add metadata
            ds = ds.sortby(ds.station_id)
            attrs1 = ds.attrs
            # filter metadata for station_ids in dataset
            logging.info("Writing metdata.")
            meta = df_inv.loc[df_inv["Climate ID"].isin(ds.station_id.values)]
            # rearrange column order to have lon, lat, elev first
            cols = meta.columns.tolist()
            for rr in [
                "Latitude (Decimal Degrees)",
                "Longitude (Decimal Degrees)",
                "Elevation (m)",
            ]:
                cols.remove(rr)
            cols1 = [
                "Latitude (Decimal Degrees)",
                "Longitude (Decimal Degrees)",
                "Elevation (m)",
            ]
            cols1.extend(cols)
            meta = meta[cols1]
            meta.index.rename("station", inplace=True)
            meta = meta.to_xarray()
            meta.sortby(meta["Climate ID"])
            meta = meta.assign({"station": ds.station.values})

            meta = meta.drop(
                ["Longitude", "Latitude"]
            )  # these values are projected x,y values Need to know prj to potentially rename
            np.testing.assert_array_equal(
                meta["Climate ID"].values, ds.station_id.values
            )
            ds = xr.merge([ds, meta])
            ds.attrs = attrs1
            del attrs1

            # TODO rename Longitude / Latitude DD
            rename = {
                "Latitude (Decimal Degrees)": "lat",
                "Longitude (Decimal Degrees)": "lon",
            }
            for i in rename.items():
                ds = ds.rename({i[0]: i[1]})

            valid_stations = list(sorted(stations_to_keep))
            valid_stations_count = len(valid_stations)

            logging.info(f"Processing stations for variable `{variable_name}`.")

            if len(station_file_codes) == 0:
                logging.error(
                    "No stations were found containing variable filename `{variable_name}`. Exiting."
                )
                return

            logging.warning(
                "Files exist for {} ECCC stations. Metadata found for {} stations. Rejecting {} stations.".format(
                    len(station_file_codes),
                    valid_stations_count,
                    len(rejected_stations),
                )
            )
            logging.warning(
                "Rejected station codes are the following: {}.".format(
                    ", ".join(rejected_stations)
                )
            )

            # Find the time dimensions for all the files
            # list_years = set()

            if hourly:
                # for i, s in enumerate(valid_stations):
                #     files = [
                #         int(Path(f).stem.split("_")[-1])
                #         for f in Path(source_files).rglob(
                #             "{}*{}*{}*.nc".format(s, variable_code, variable_name)
                #         )
                #     ]
                #     list_years.update(files)

                # list_years = [int(Path(f).stem.split("_")[-1]) for f in list_files_to_combine]
                year_start = ds.time.dt.year.min().values
                year_end = ds.time.dt.year.max().values
                # year_start, year_end = min(list_years), max(list_years)

                # Calculate the dimensions of the output NetCDF
                time_index = pd.date_range(
                    start="{}-01-01".format(year_start),
                    end="{}-01-01".format(year_end + 1),
                    freq="H",
                )[:-1]

            else:
                time_index = pd.date_range(
                    start="1840-01-01", end="{}-01-01".format(dt.today().year), freq="D"
                )

            logging.info("Preparing the NetCDF.")
            logging.info(
                "Number of ECCC stations: {}, time steps: {}.".format(
                    valid_stations_count, time_index.size
                )
            )

            dsOut = xr.Dataset(
                coords={
                    "time": time_index,
                    "station": ds.station,
                    "station_id": ds.station_id,
                },
                attrs=ds.attrs,
            )

            for vv in ds.data_vars:
                dsOut[vv] = ds[
                    vv
                ]  # assign data variables to output dataset ... will align with time coords

            output_folder = output_folder.joinpath("merged")
            output_folder.mkdir(parents=True, exist_ok=True)

            file_out = Path(output_folder).joinpath(
                "{}_eccc_{}".format(variable_name, "hourly" if hourly else "daily",)
            )

            if mf_dataset_freq is not None:
                _, datasets = zip(
                    *dsOut.resample(time=mf_dataset_freq)
                )  # output mf_dataseset using resampling frequency
            else:
                datasets = [dsOut]

            paths = [
                f'{file_out}_{dd.time.dt.year.min().values}-{dd.time.dt.year.max().values}_created{dt.now().strftime("%Y%m%d")}.nc'
                for dd in datasets
            ]

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in datasets[0].data_vars}

            with ProgressBar():
                for ii in zip(datasets, paths):
                    dd, path = ii
                    dd.to_netcdf(
                        path, engine="h5netcdf", format="NETCDF4", encoding=encoding
                    )
        else:
            logging.info("No files found for variable `{}`.".format(variable_name))

    logging.warning(
        "Process completed in {:.2f} seconds".format(time.time() - func_time)
    )

    #     if file_out.exists():
    #         file_out.unlink()
    #
    #     ds = netCDF4.Dataset(file_out, "w", format="NETCDF4")
    #     ds.createDimension("time", None)
    #     ds.createDimension("station", valid_stations_count)
    #
    #     # creation de la variable
    #     least_significant_digit = info["least_significant_digit"]
    #     nc_var = ds.createVariable(
    #         variable_name,
    #         datatype="f4",
    #         dimensions=("time", "station"),
    #         zlib=True,
    #         least_significant_digit=least_significant_digit,
    #         chunksizes=(100000, 1),
    #     )
    #     nc_var.units = info["nc_units"]
    #
    #     nc_flag = None
    #     if include_flags:
    #         # Create variable for the flags
    #         nc_flag = ds.createVariable(
    #             "flag",
    #             datatype="S1",
    #             dimensions=("time", "station"),
    #             zlib=True,
    #             chunksizes=(100000, 1),
    #         )
    #         nc_flag.long_name = "data flag"
    #         nc_flag.note = "See ECCC technical documentation for details"
    #
    #     # Create the time variable
    #     nc_time = ds.createVariable("time", "f8", dimensions="time")
    #     tunits = time_index[0].strftime("days since %Y-%m-%d %H:%M:%S")
    #     nc_time.units = tunits
    #     nc_time.calendar = "proleptic_gregorian"
    #     time_num = cftime.date2num(time_index.to_pydatetime(), tunits)
    #     nc_time[:] = time_num
    #     ds.sync()
    #
    #     # Write out the appropriate ECCC metadata
    #     nc_name = ds.createVariable("name", "str", dimensions="station")
    #     nc_name.long_name = "Name"
    #     nc_lon = ds.createVariable("lon", "f8", dimensions="station")
    #     nc_lon.units = "degrees"
    #     nc_lon.long_name = "Longitude (Decimal Degrees)"
    #     nc_lon.standard_name = "longitude"
    #     nc_lat = ds.createVariable("lat", "f8", dimensions="station")
    #     nc_lat.units = "degrees"
    #     nc_lat.long_name = "Latitude (Decimal Degrees)"
    #     nc_lat.standard_name = "latitude"
    #     nc_province = ds.createVariable("province", "str", dimensions="station")
    #     nc_province.long_name = "Province"
    #     nc_elevation = ds.createVariable("elevation", "f", dimensions="station")
    #     nc_elevation.units = "m"
    #     nc_elevation.long_name = "Elevation"
    #     nc_elevation.standard_name = "height"
    #     nc_climate_id = ds.createVariable("climate_id", "str", dimensions="station")
    #     nc_climate_id.long_name = "Climate ID"
    #     nc_station_id = ds.createVariable("station_id", "str", dimensions="station")
    #     nc_station_id.long_name = "Station ID"
    #     nc_wmo_id = ds.createVariable("wmo_id", "str", dimensions="station")
    #     nc_wmo_id.long_name = "WMO ID"
    #     nc_tc_id = ds.createVariable("tc_id", "str", dimensions="station")
    #     nc_tc_id.long_name = "TC ID"
    #     ds.sync()
    #
    #     # Use a Pandas DataFrame to align the data
    #     # df_tot = pd.DataFrame(index=time_index)
    #     # df_tot = xr.Dataset(coords=dict(time=time_index))
    #
    #     # Grab ECCC station metadata for entry
    #     for iter_station, station in enumerate(valid_stations):
    #         t00 = time.time()
    #         logging.info(
    #             "Collecting data from station {}/{}: Station ID {}.".format(
    #                 iter_station + 1, valid_stations_count, station
    #             )
    #         )
    #         df_stat = df_inv.loc[df_inv["Climate ID"] == station]
    #
    #         # Only take one metadata entry per station
    #         if df_stat.shape[0] != 1:
    #             logging.warning(
    #                 "Duplicate metadata found for station {}. Skipping.".format(station)
    #             )
    #             continue
    #         df_stat = df_stat.iloc[0]
    #
    #         # Open files with special handling for multi-year data
    #         list_station_files = [
    #             f
    #             for f in source_files.rglob(
    #                 "{}*{}*{}*.nc".format(station, variable_code, variable_name)
    #             )
    #         ]
    #         # single_year_files = []
    #         # multi_year_files = []
    #         # for file in list_station_files:
    #         #     a1 = file.name.split("_")[-2]
    #         #     a2 = file.name.split("_")[-1][:4]
    #         #     if a1 != a2:
    #         #         multi_year_files.append(file)
    #         #     else:
    #         #         single_year_files.append(file)
    #
    #         # Annual data treatment
    #         if len(list_station_files) > 0:
    #             try:
    #                 ds_single = xr.open_mfdataset(
    #                     list_station_files, combine="by_coords"
    #                 )
    #             except ValueError:
    #                 logging.exception(
    #                     "Unable to read by coordinates. Skipping station: {}".format(
    #                         station
    #                     )
    #                 )
    #                 continue
    #         else:
    #             continue
    #             #
    #             # df_tot[variable_name] = ds_single[variable_name].to_dataframe()
    #             #
    #             # if include_flags:
    #             #     df_tot["flag"] = ds_single["flag"].to_dataframe()
    #
    #         # Add multi-annual data
    #         # for fichier in multi_year_files:
    #         #     logging.info(
    #         #         "Special handling for multi-year data file: {}".format(fichier)
    #         #     )
    #         #     # Read data into a dataframe and keep valid entries only
    #         #     df = xr.open_dataset(fichier).to_dataframe()
    #         #     df_valid = df.loc[~df[variable_name].isnull()]
    #         #
    #         #     # Traitement special pour les donnees en double
    #         #     #
    #         #     # Dans les fichiers multi-annuels, il arrive qu'il y ait plusieurs
    #         #     # entrees pour les memes dates/eccc
    #         #     #
    #         #     # *** CA NE DEVRAIT PAS ARRIVER ***
    #         #     #
    #         #     if not df_valid.index.is_unique:
    #         #         logging.warning("Donnees en double ... on applique une patch ...")
    #         #         # Pour les temps ou on a plus d'une donnee valide, on garde la derniere
    #         #         df_valid = df_valid[
    #         #             not df_valid.index.duplicated(keep=double_handling)
    #         #         ]
    #         #         assert df_valid.index.is_unique
    #         #
    #         #     if df_valid.index.size != 0:
    #         #         df_tot.loc[df_valid.index, variable_name] = df_valid[variable_name]
    #         #         if include_flags:
    #         #             df_tot.loc[df_valid.index, "flag"] = df_valid["flag"]
    #
    #         # Dump the data into the output container
    #         # nc_var[:, iter_station] = df_tot[variable_name].values
    #         logging.info("Or maybe it's over here?")
    #         nc_var[:, iter_station] = ds_single[variable_name].values
    #         logging.info("It could be this?")
    #
    #         if include_flags:
    #             # nc_flag[:, iter_station] = df_tot["flag"].values
    #             nc_flag[:, iter_station] = ds_single["flag"].values
    #
    #         logging.info(
    #             "Added {:} data entries.".format(int(ds_single[variable_name].count()))
    #         )
    #
    #         # # Remove the df_tot values
    #         # df_tot.drop(columns=[variable_name])
    #         # if include_flags:
    #         #     df_tot.drop(columns=["flag"])
    #
    #         # Add the metadata to the file
    #         nc_name[iter_station] = df_stat["Name"]
    #         nc_province[iter_station] = df_stat["Province"]
    #         nc_lat[iter_station] = df_stat["Latitude (Decimal Degrees)"]
    #         nc_lon[iter_station] = df_stat["Longitude (Decimal Degrees)"]
    #         nc_elevation[iter_station] = df_stat["Elevation (m)"]
    #         nc_climate_id[iter_station] = str(df_stat["Climate ID"])
    #         nc_station_id[iter_station] = str(df_stat["Station ID"])
    #         nc_wmo_id[iter_station] = str(df_stat["WMO ID"])
    #         nc_tc_id[iter_station] = str(df_stat["TC ID"])
    #         ds.sync()
    #
    #         logging.info(
    #             "Transaction duration: {:.2f} seconds.".format(time.time() - t00)
    #         )
    #
    #     ds.Conventions = "CF-1.7"
    #
    #     ds.title = "Environment and Climate Change Canada (ECCC) weather eccc"
    #     ds.history = "{}: Merged from multiple individual station files to n-dimensional array.".format(
    #         dt.now().strftime("%Y-%m-%d %X")
    #     )
    #     ds.version = "v{}".format(dt.now().strftime("%Y.%m"))
    #     ds.institution = "Environment and Climate Change Canada (ECCC)"
    #     ds.source = "Weather Station data <ec.services.climatiques-climate.services.ec@canada.ca>"
    #     ds.references = "https://climate.weather.gc.ca/doc/Technical_Documentation.pdf"
    #     ds.comment = "Acquired on demand from data specialists at ECCC Climate Services / Services Climatiques"
    #     ds.redistribution = "Redistribution policy unknown. For internal use only."
    #
    #     ds.close()
    # logging.warning(
    #     "Process completed in {:.2f} seconds".format(time.time() - func_time)
    # )
