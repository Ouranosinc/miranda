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
import itertools
import logging
import tempfile
import time
from calendar import monthrange
from datetime import datetime as dt
from logging import config
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar

from miranda.scripting import LOGGING_CONFIG

from ._utils import cf_daily_metadata, cf_hourly_metadata

config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "aggregate_stations",
    "convert_hourly_flat_files",
    "convert_daily_flat_files",
    "merge_converted_variables",
]


def convert_hourly_flat_files(
    source_files: Union[str, Path],
    output_folder: Union[str, Path, List[Union[str, int]]],
    variables: Union[str, int, List[Union[str, int]]],
    missing_value: int = -9999,
) -> None:
    """

    Parameters
    ----------
    source_files : str or Path
    output_folder : str or Path
    variables : str or List[str]
    missing_value : int

    Returns
    -------
    None
    """
    func_time = time.time()

    if isinstance(variables, (str, int)):
        variables = [variables]

    for variable_code in variables:
        info = cf_hourly_metadata(variable_code)
        variable_code = str(variable_code).zfill(3)
        variable_name = info["standard_name"]
        variable_file_name = info["nc_name"]

        # Preparing the data extraction
        col_names = "code year month day code_var ".split()
        for i in range(1, 25):
            col_names.append(f"D{i:0n}")
            col_names.append(f"F{i:0n}")

        rep_nc = Path(output_folder).joinpath(variable_file_name)
        rep_nc.mkdir(parents=True, exist_ok=True)

        # Loop on the files
        list_files = list()
        if isinstance(source_files, list) or Path(source_files).is_file():
            list_files.append(source_files)
        elif 262 < int(variable_code) <= 280:
            list_files.extend(
                [f for f in Path(source_files).rglob("HLY*RCS*") if f.is_file()]
            )
        else:
            list_files.extend(
                [f for f in Path(source_files).rglob("HLY*") if f.is_file()]
            )

        errored_files = list()
        for fichier in list_files:
            logging.info(f"Processing file: {fichier}.")

            # Create a dataframe from the files
            try:
                df = pd.read_fwf(
                    fichier,
                    widths=[7, 4, 2, 2, 3] + [6, 1] * 24,
                    names=col_names,
                    dtype={"year": int, "month": int, "day": int, "code_var": str},
                )
            except FileNotFoundError:
                logging.error(f"File {fichier} was not found.")
                errored_files.append(fichier)
                continue

            except (UnicodeDecodeError, Exception):
                logging.error(
                    f"File {fichier} was unable to be read. This is probably an issue with the file."
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
                dfd = df_var.loc[:, [f"D{i:0n}" for i in range(1, 25)]]
                dff = df_var.loc[:, [f"F{i:0n}" for i in range(1, 25)]]

                # Remove the "NaN" flag
                dff = dff.fillna("")

                # Use the flag to mask the values
                try:
                    val = np.asfarray(dfd.values)
                except ValueError as e:
                    logging.error(f"{e} raised from {dfd}, continuing...")
                    continue
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
                ds.attrs["version"] = f"v{dt.now().strftime('%Y.%m')}"
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

    logging.warning(f"Process completed in {time.time() - func_time:.2f} seconds")


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
        info = cf_daily_metadata(variable_code)
        variable_code = str(variable_code).zfill(3)
        nc_name = info["nc_name"]

        # Prepare the data extraction
        titre_colonnes = "code year month code_var".split()
        for i in range(1, 32):
            titre_colonnes.append(f"D{i:0n}")
            titre_colonnes.append(f"F{i:0n}")

        # Create the output directory
        rep_nc = Path(output_folder).joinpath(nc_name)
        rep_nc.mkdir(parents=True, exist_ok=True)

        # Loop on the files
        list_files = list()
        if isinstance(source_files, list) or Path(source_files).is_file():
            list_files.append(source_files)
        else:
            list_files.extend(
                [f for f in Path(source_files).rglob("*DLY*") if f.is_file()]
            )

        errored_files = list()
        for fichier in list_files:
            logging.info(f"Processing file: {fichier}.")

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
                logging.info(f"Converting {nc_name} for station code: {code}")

                # Dump the values into a DataFrame
                df_var = df_code[df_code["code_var"] == variable_code].copy()

                # Apply the mask according to the NaN value
                df_var = df_var.replace(missing_value, np.nan)

                # Decoding the values and flags
                dfd = df_var.loc[:, [f"D{i:0n}" for i in range(1, 32)]]
                dff = df_var.loc[:, [f"F{i:0n}" for i in range(1, 32)]]

                # Remove the "NaN" flag
                dff = dff.fillna("")

                try:
                    # Use the flag to mask the values
                    val = np.asfarray(dfd.values)
                    flag = dff.values
                    mask = np.isin(flag, info["missing_flags"])
                    val[mask] = np.nan
                except ValueError:
                    continue

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
                        c=code,
                        vc=variable_code,
                        v=nc_name,
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

    logging.warning(f"Process completed in {time.time() - func_time:.2f} seconds")


def aggregate_stations(
    source_files: Optional[Union[str, Path]] = None,
    output_folder: Optional[Union[str, Path]] = None,
    station_metadata: Union[str, Path] = None,
    time_step: str = "h",
    variables: Optional[Union[str, int, List[Union[str, int]]]] = None,
    include_flags: bool = True,
    groups: int = 5,
    mf_dataset_freq: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
) -> None:
    """

    Parameters
    ----------
    source_files: Union[str, Path]
    output_folder: Union[str, Path]
    variables: Optional[Union[str, int, List[Union[str, int]]]]
    time_step: str
    station_metadata: Union[str, Path]
    include_flags: bool
    groups: int
      The number of file groupings used for converting to multi-file Datasets.
    mf_dataset_freq: Optional[str]
      Resampling frequency for creating output multi-file Datasets. E.g. 'YS': 1 year per file, '5YS': 5 years per file.
    temp_directory: Optional[Union[str, Path]]
      Use another temporary directory location in case default location is not spacious enough.

    Returns
    -------
    None
    """
    func_time = time.time()

    if not station_metadata:
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
            info = cf_hourly_metadata(variable_code)
        else:
            info = cf_daily_metadata(variable_code)
        variable_name = info["nc_name"]
        logging.info(f"Merging `{variable_name}` using `{time_step}` time step.")

        # Find the ECCC stations where we have available metadata
        df_inv = pd.read_csv(str(station_metadata), header=3)
        station_inventory = list(df_inv["Climate ID"].values)

        # Only perform aggregation on available data with corresponding metadata
        logging.info("Performing glob and sort.")
        nclist = sorted(list(source_files.joinpath(variable_name).rglob("*.nc")))

        ds = None
        if nclist != list():
            nclists = np.array_split(nclist, groups)

            with tempfile.TemporaryDirectory(
                prefix="eccc", dir=temp_directory
            ) as temp_dir:
                combinations = [(ii, nc, temp_dir) for ii, nc in enumerate(nclists)]

                # TODO memory use seems ok here .. could try using Pool() to increase performance
                for combo in combinations:
                    ii, nc, temp_dir = combo
                    _tmp_nc(ii, nc, temp_dir, groups)

                ds = xr.open_mfdataset(
                    sorted(list(Path(temp_dir).glob("*.nc"))),
                    combine="nested",
                    concat_dim="station",
                    chunks=dict(time=365),
                )

                # dask gives warnings about export 'object' data types
                ds["station_id"] = ds["station_id"].astype(str)
        if ds:
            station_file_codes = [x.name.split("_")[0] for x in nclist]
            rejected_stations = set(station_file_codes).difference(
                set(station_inventory)
            )

            logging.info(f"{len(rejected_stations)} rejected due to missing metadata.")
            r_all = np.zeros(ds.station_id.shape) == 0

            for r in rejected_stations:
                r_all[ds.station_id == r] = False
            ds = ds.isel(station=r_all)
            if not include_flags:
                drop_vars = [vv for vv in ds.data_vars if "flag" in vv]
                ds = ds.drop_vars(drop_vars)

            # Ensure data is in order to add metadata
            ds = ds.sortby(ds.station_id)

            attrs1 = ds.attrs
            # filter metadata for station_ids in dataset
            logging.info("Writing metdata.")

            meta = df_inv.loc[df_inv["Climate ID"].isin(ds.station_id.values)]
            # Rearrange column order to have lon, lat, elev first
            cols = meta.columns.tolist()
            cols1 = [
                "Latitude (Decimal Degrees)",
                "Longitude (Decimal Degrees)",
                "Elevation (m)",
            ]
            for rr in cols1:
                cols.remove(rr)
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

            # TODO rename Longitude / Latitude DD
            rename = {
                "Latitude (Decimal Degrees)": "lat",
                "Longitude (Decimal Degrees)": "lon",
            }
            for i in rename.items():
                ds = ds.rename({i[0]: i[1]})

            valid_stations = list(sorted(ds.station_id.values))
            valid_stations_count = len(valid_stations)

            logging.info(f"Processing stations for variable `{variable_name}`.")

            if len(station_file_codes) == 0:
                logging.error(
                    f"No stations were found containing variable filename `{variable_name}`. Exiting."
                )
                return

            logging.warning(
                "Files exist for {} ECCC stations. Metadata found for {} stations. Rejecting {} stations.".format(
                    len(station_file_codes),
                    valid_stations_count,
                    len(rejected_stations),
                )
            )
            if rejected_stations:
                logging.warning(
                    f"Rejected station codes are the following: {', '.join(rejected_stations)}."
                )

            logging.info("Preparing the NetCDF time period.")
            # Create the time period timestamps
            year_start = ds.time.dt.year.min().values
            year_end = ds.time.dt.year.max().values

            # Calculate the time index dimensions of the output NetCDF
            time_index = pd.date_range(
                start=f"{year_start}-01-01",
                end=f"{year_end + 1}-01-01",
                freq="H" if hourly else "D",
            )[:-1]

            logging.info(
                "Number of ECCC stations: {}, time steps: {}.".format(
                    valid_stations_count, time_index.size
                )
            )

            ds_out = xr.Dataset(
                coords={
                    "time": time_index,
                    "station": ds.station,
                    "station_id": ds.station_id,
                },
                attrs=ds.attrs,
            )

            for vv in ds.data_vars:
                ds_out[vv] = ds[
                    vv
                ]  # assign data variables to output dataset ... will align with time coords

            output_folder.mkdir(parents=True, exist_ok=True)

            file_out = Path(output_folder).joinpath(
                "{}_eccc_{}".format(
                    variable_name,
                    "hourly" if hourly else "daily",
                )
            )

            if mf_dataset_freq is not None:
                _, datasets = zip(
                    *ds_out.resample(time=mf_dataset_freq)
                )  # output mf_dataset using resampling frequency
            else:
                datasets = [ds_out]

            paths = [
                f"{file_out}_{dd.time.dt.year.min().values}-{dd.time.dt.year.max().values}_"
                f'created{dt.now().strftime("%Y%m%d")}.nc'
                for dd in datasets
            ]

            comp = dict(zlib=True, complevel=5)

            with ProgressBar():
                for dataset, path in zip(datasets, paths):
                    encoding = {var: comp for var in ds_out.data_vars}
                    dataset.to_netcdf(
                        path,
                        engine="h5netcdf",
                        format="NETCDF4",
                        encoding=encoding,
                    )
                    dataset.close()
                    del dataset
            ds.close()
            ds_out.close()

        else:
            logging.info("No files found for variable: `%s`." % variable_name)

    runtime = f"Process completed in {time.time() - func_time:.2f} seconds"
    logging.warning(runtime)


def _tmp_nc(
    ii: int,
    nc: Union[str, Path],
    tempdir: Union[str, Path],
    batches: Optional[int] = None,
) -> None:
    if batches is None:
        batches = "X"
    logging.info(f"Processing batch of files {ii + 1} of {batches}")
    station_file_codes = [x.name.split("_")[0] for x in nc]

    ds = xr.open_mfdataset(nc, combine="nested", concat_dim="station")
    ds = ds.assign_coords(
        station_id=xr.DataArray(station_file_codes, dims="station").astype(str)
    )
    if "flag" in ds.data_vars:
        ds1 = ds.drop_vars("flag").copy(deep=True)
        ds1["flag"] = ds.flag.astype(str)
        ds = ds1

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}

    with ProgressBar():
        ds.load().to_netcdf(
            Path(tempdir).joinpath(f"{str(ii).zfill(3)}.nc"),
            engine="h5netcdf",
            format="NETCDF4",
            encoding=encoding,
        )
        del ds


def merge_converted_variables(
    source: Union[str, Path],
    destination: Union[str, Path],
    variables: Optional[Union[str, int, List[Union[str, int]]]] = None,
) -> None:
    """

    Parameters
    ----------
    source : Union[str, Path]
    destination : Union[str, Path]
    variables : Optional[Union[str, int, List[Union[str, int]]]]

    Returns
    -------

    """

    def _combine_years(args: Tuple[str, Union[str, Path], Union[str, Path]]) -> None:
        varia, input_folder, output_folder = args

        ncfiles = sorted(list(input_folder.glob("*.nc")))
        logging.info(
            f"Found {len(ncfiles)} files for station code {input_folder.name}."
        )
        logging.info(f"Opening: {ncfiles}")

        ds = xr.open_mfdataset(
            ncfiles, parallel=False, combine="by_coords", concat_dim={"time"}
        )

        outfile = output_folder.joinpath(
            f'{ncfiles[0].name.split(f"_{varia}_")[0]}_{varia}_'
            f"{ds.time.dt.year.min().values}-{ds.time.dt.year.max().values}.nc"
        )
        if not outfile.exists():
            logging.info(f"Merging to {outfile.name}")
            comp = dict(zlib=True, complevel=5)
            encoding = {data_var: comp for data_var in ds.data_vars}
            encoding["time"] = {"dtype": "single"}
            with ProgressBar():
                ds.to_netcdf(outfile, encoding=encoding)
        else:
            logging.info(f"Files exist for {outfile.name}. Continuing...")

    if isinstance(source, str):
        source = Path(source)
    if isinstance(destination, str):
        destination = Path(destination)

    selected_variables = list()
    if variables is not None:
        if not isinstance(variables, list):
            variables = [variables]
        for var in variables:
            try:
                selected_variables.append(cf_hourly_metadata(var))
            except KeyError:
                selected_variables.append(cf_hourly_metadata(var))

    variables_found = [x.name for x in source.iterdir() if x.is_dir()]
    if selected_variables:
        variables_found = [
            x
            for x in variables_found
            if x in [item["nc_name"] for item in selected_variables]
        ]

    for variable in variables_found:
        logging.info(f"Merging files found for variable: `{variable}`.")
        station_dirs = [x for x in source.joinpath(variable).iterdir() if x.is_dir()]
        logging.info(f"Number of stations found: {len(station_dirs)}.")
        outrep = destination.joinpath(variable)

        Path(outrep).mkdir(parents=True, exist_ok=True)
        combs = list(itertools.product(*[[variable], station_dirs, [outrep]]))
        for c in combs:
            try:
                _combine_years(c)
            except ValueError as e:
                logging.error(
                    f"`{e}` encountered for station `{c[1].name}`. Continuing..."
                )
