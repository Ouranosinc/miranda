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
import contextlib
import functools
import logging
import multiprocessing as mp
import os
import re
import sys
import tempfile
import time
from calendar import monthrange
from datetime import datetime as dt
from logging import config
from pathlib import Path
from typing import List, Optional, Set, Union
from urllib.error import HTTPError

import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar
from xclim.core.units import convert_units_to

from miranda.scripting import LOGGING_CONFIG
from miranda.storage import file_size, report_file_size
from miranda.units import GiB, MiB
from miranda.utils import generic_extract_archive

from ._utils import cf_station_metadata

config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "aggregate_stations",
    "convert_flat_files",
    "merge_converted_variables",
]

TABLE_DATE = dt.now().strftime("%d %B %Y")


def load_station_metadata(meta: Union[str, Path]) -> xr.Dataset:
    if meta:
        df_inv = pd.read_csv(meta, header=0)
    else:
        try:
            import geopandas as gpd

            station_metadata_url = "https://api.weather.gc.ca/collections/climate-stations/items?f=json&limit=15000000"
            df_inv = gpd.read_file(station_metadata_url)
        except HTTPError as err:
            raise RuntimeError(
                f"Station metadata table unable to be fetched. Considering downloading directly: {err}"
            )
    df_inv["LONGITUDE"] = df_inv.geometry.x
    df_inv["LATITUDE"] = df_inv.geometry.y
    df_inv["ELEVATION"] = df_inv.ELEVATION.astype(float)
    df_inv["CLIMATE_IDENTIFIER"] = df_inv["CLIMATE_IDENTIFIER"].astype(str)

    df_inv = df_inv.drop(["geometry"], axis=1)
    return df_inv.to_xarray()


def _convert_station_file(
    fichier: Path,
    output_path: Path,
    errored_files: List[Path],
    mode: str,
    add_offset: float,
    column_dtypes: List[str],
    column_names: List[str],
    long_name: str,
    missing_flags: Set[str],
    missing_values: Set[str],
    nc_name: str,
    raw_units: str,
    units: str,
    scale_factor: float,
    standard_name: str,
    variable_code: str,
    **kwargs,
):
    if mode.lower() in ["h", "hour", "hourly"]:
        num_observations = 24
        column_widths = [7, 4, 2, 2, 3] + [6, 1] * num_observations
    elif mode.lower() in ["d", "day", "daily"]:
        num_observations = 31
        column_widths = [7, 4, 2, 3] + [6, 1] * num_observations
    else:
        raise NotImplementedError("`mode` must be 'h'/'hourly or 'd'/'daily'.")

    if not missing_values:
        missing_values = {-9999, "#####"}

    with tempfile.TemporaryDirectory() as temp_folder:
        if fichier.suffix in [".gz", ".tar", ".zip", ".7z"]:
            data_files = generic_extract_archive(fichier, output_dir=temp_folder)
        else:
            data_files = [fichier]
        logging.info(f"Processing file: {fichier}.")

        size_limit = 1 * GiB

        for data in data_files:
            if file_size(data) > size_limit and "dask" in sys.modules:
                logging.info(
                    f"File exceeds {report_file_size(size_limit)} - Using dask.dataframes."
                )
                pandas_reader = dd
                using_dask_array = True
                chunks = dict(blocksize=200 * MiB)
                client = ProgressBar
            else:
                logging.info(
                    f"File below {report_file_size(size_limit)} - Using pandas.dataframes."
                )
                pandas_reader = pd
                chunks = dict()
                using_dask_array = False
                client = contextlib.nullcontext

            with client():
                # Create a dataframe from the files
                try:
                    df = pandas_reader.read_fwf(
                        data,
                        widths=column_widths,
                        names=column_names,
                        dtype={
                            name: data_type
                            for name, data_type in zip(column_names, column_dtypes)
                        },
                        assume_missing=True,
                        **chunks,
                    )
                except FileNotFoundError:
                    logging.error(f"File {data} was not found.")
                    errored_files.append(data)
                    return

                except (UnicodeDecodeError, Exception) as e:
                    logging.error(
                        f"File {data.name} was unable to be read. "
                        f"This is probably an issue with the file: {e}"
                    )
                    errored_files.append(data)
                    return

                # Loop through the station codes
                station_codes = df["code"].unique()
                for code in station_codes:
                    df_code = df[df["code"] == code]

                    # Abort if the variable is not found
                    if using_dask_array:
                        has_variable_codes = (
                            (df_code["code_var"] == variable_code).compute()
                        ).any()
                    else:
                        has_variable_codes = (
                            df_code["code_var"] == variable_code
                        ).any()
                    if not has_variable_codes:
                        logging.info(
                            f"Variable `{nc_name}` not found for station code: {code} in file {data}. Continuing..."
                        )
                        continue

                    # Perform the data treatment
                    logging.info(f"Converting `{nc_name}` for station code: {code}")

                    # Dump the data into a DataFrame
                    df_var = df_code[df_code["code_var"] == variable_code].copy()

                    # Mask the data according to the missing values flag
                    df_var = df_var.replace(missing_values, np.nan)

                    # Decode the values and flags
                    dfd = df_var.loc[
                        :, [f"D{i:0n}" for i in range(1, num_observations + 1)]
                    ]
                    dff = df_var.loc[
                        :, [f"F{i:0n}" for i in range(1, num_observations + 1)]
                    ]

                    # Remove the "NaN" flag
                    dff = dff.fillna("")

                    # Use the flag to mask the values
                    try:
                        val = np.asarray(dfd.values, float)
                    except ValueError as e:
                        logging.error(f"{e} raised from {dfd}, continuing...")
                        continue
                    try:
                        flag = np.asarray(dff.values, str)
                    except ValueError as e:
                        logging.error(f"{e} raised from {dff}, continuing...")
                        continue
                    mask = np.isin(flag, missing_flags)
                    val[mask] = np.nan

                    # Treat according to units conversions
                    val = val * scale_factor + add_offset

                    # Create the DataArray
                    date_summations = dict(time=list())
                    if mode == "hourly":
                        for index, row in df_var.iterrows():
                            period = pd.Period(
                                year=row.year, month=row.month, day=row.day, freq="D"
                            )
                            dates = pd.Series(
                                pd.date_range(
                                    start=period.start_time,
                                    end=period.end_time,
                                    freq="H",
                                )
                            )
                            date_summations["time"].extend(dates)
                        written_values = val.flatten()
                        written_flags = flag.flatten()
                    elif mode == "daily":
                        value_days = list()
                        flag_days = list()
                        for i, (index, row) in enumerate(df_var.iterrows()):
                            period = pd.Period(year=row.year, month=row.month, freq="M")
                            dates = pd.Series(
                                pd.date_range(
                                    start=period.start_time,
                                    end=period.end_time,
                                    freq="D",
                                )
                            )
                            date_summations["time"].extend(dates)

                            value_days.extend(
                                val[i][
                                    range(monthrange(int(row.year), int(row.month))[1])
                                ]
                            )
                            flag_days.extend(
                                flag[i][
                                    range(monthrange(int(row.year), int(row.month))[1])
                                ]
                            )
                        written_values = value_days
                        written_flags = flag_days

                    ds = xr.Dataset()
                    da_val = xr.DataArray(
                        written_values, coords=date_summations, dims=["time"]
                    )

                    if raw_units != units:
                        da_val.attrs["units"] = raw_units
                        da_val = convert_units_to(da_val, units)
                    else:
                        da_val.attrs["units"] = units

                    da_val = da_val.rename(nc_name)
                    variable_attributes = dict(
                        variable_code=variable_code,
                        standard_name=standard_name,
                        long_name=long_name,
                    )
                    if "original_units" in kwargs:
                        variable_attributes["original_units"] = kwargs["original_units"]
                    da_val.attrs.update(variable_attributes)

                    da_flag = xr.DataArray(
                        written_flags, coords=date_summations, dims=["time"]
                    )
                    da_flag = da_flag.rename("flag")
                    flag_attributes = dict(
                        long_name="data flag",
                        note="See ECCC technical documentation for details",
                    )
                    da_flag.attrs.update(flag_attributes)

                    ds[nc_name] = da_val
                    ds["flag"] = da_flag

                    # save the file in NetCDF format
                    start_year = ds.time.dt.year.values[0]
                    end_year = ds.time.dt.year.values[-1]

                    station_folder = output_path.joinpath(str(code))
                    station_folder.mkdir(parents=True, exist_ok=True)

                    f_nc = (
                        f"{code}_{variable_code}_{nc_name}_"
                        f"{start_year if start_year == end_year else '_'.join([str(start_year), str(end_year)])}.nc"
                    )

                    if station_folder.joinpath(f_nc).exists():
                        logging.warning(f"File `{f_nc}` already exists. Continuing...")

                    history = (
                        f"{dt.now().strftime('%Y-%m-%d %X')} converted from flat station file "
                        f"(`{fichier.name}`) to n-dimensional array."
                    )

                    # TODO: This info should eventually be sourced from a JSON definition
                    global_attrs = dict(
                        Conventions="CF-1.8",
                        comment="Acquired on demand from data specialists at "
                        "ECCC Climate Services / Services Climatiques.",
                        contact="John Richard",
                        contact_email="climatcentre-climatecentral@ec.gc.ca",
                        domain="CAN",
                    )
                    if mode == "hourly":
                        global_attrs.update(dict(frequency="1hr"))
                    elif mode == "daily":
                        global_attrs.update(dict(frequency="day"))
                    global_attrs.update(
                        dict(
                            history=history,
                            institution="ECCC",
                            member=code,
                            processing_level="raw",
                            redistribution="Redistribution policy unknown. For internal use only.",
                            references="https://climate.weather.gc.ca/doc/Technical_Documentation.pdf",
                            source="historical-station-records",
                            table_date=TABLE_DATE,
                            title="Environment and Climate Change Canada (ECCC) weather eccc",
                            type="station-obs",
                            variable=str(nc_name),
                            version=f"v{dt.now().strftime('%Y.%m.%V')}",  # Year.Month.Week
                        )
                    )
                    ds.attrs.update(global_attrs)

                    logging.info(f"Exporting to: {station_folder.joinpath(f_nc)}")
                    ds.to_netcdf(station_folder.joinpath(f_nc))
                    del ds
                    del val
                    del mask
                    del flag
                    del da_val
                    del da_flag
                    del dfd
                    del dff
                    del written_values
                    del written_flags
                    del date_summations

                del df

        if os.listdir(temp_folder):
            for temporary_file in Path(temp_folder).glob("*"):
                if temporary_file in data_files:
                    temporary_file.unlink()


def convert_flat_files(
    source_files: Union[str, Path],
    output_folder: Union[str, Path, List[Union[str, int]]],
    variables: Union[str, int, List[Union[str, int]]],
    mode: str = "hourly",
    n_workers: int = 4,
) -> None:
    """

    Parameters
    ----------
    source_files: str or Path
    output_folder: str or Path
    variables: str or List[str]
    mode: {"hourly", "daily"}
    n_workers: int

    Returns
    -------
    None
    """
    func_time = time.time()

    if mode.lower() in ["h", "hour", "hourly"]:
        num_observations = 24
        column_names = ["code", "year", "month", "day", "code_var"]
        column_dtypes = [str, float, float, float, str]
    elif mode.lower() in ["d", "day", "daily"]:
        num_observations = 31
        column_names = ["code", "year", "month", "code_var"]
        column_dtypes = [str, float, float, str]
    else:
        raise NotImplementedError("`mode` must be 'h'/'hourly or 'd'/'daily'.")

    # Preparing the data column headers
    for i in range(1, num_observations + 1):
        data_entry, flag_entry = f"D{i:0n}", f"F{i:0n}"
        column_names.append(data_entry)
        column_names.append(flag_entry)
        column_dtypes.extend([str, str])

    if isinstance(variables, (str, int)):
        variables = [variables]

    for variable_code in variables:
        variable_code = str(variable_code).zfill(3)
        metadata = cf_station_metadata(variable_code)
        nc_name = metadata["nc_name"]

        rep_nc = Path(output_folder).joinpath(nc_name)
        rep_nc.mkdir(parents=True, exist_ok=True)

        # Loop on the files
        logging.info(
            f"Collecting files for variable '{metadata['standard_name']}' "
            f"(filenames containing '{metadata['_table_name']}')."
        )
        list_files = list()
        if isinstance(source_files, list) or Path(source_files).is_file():
            list_files.append(source_files)
        else:
            glob_patterns = [g for g in metadata["_table_name"]]
            for pattern in glob_patterns:
                list_files.extend(
                    [f for f in Path(source_files).rglob(f"{pattern}*") if f.is_file()]
                )
        manager = mp.Manager()
        errored_files = manager.list()
        converter_func = functools.partial(
            _convert_station_file,
            output_path=rep_nc,
            errored_files=errored_files,
            mode=mode,
            variable_code=variable_code,
            column_names=column_names,
            column_dtypes=column_dtypes,
            **metadata,
        )
        with mp.Pool(processes=n_workers) as pool:
            pool.map(converter_func, list_files)
            pool.close()
            pool.join()

        if errored_files:
            logging.warning(
                "Some files failed to be properly parsed:\n", ", ".join(errored_files)
            )

    logging.warning(f"Process completed in {time.time() - func_time:.2f} seconds")


def aggregate_stations(
    source_files: Optional[Union[str, Path]] = None,
    output_folder: Optional[Union[str, Path]] = None,
    time_step: str = None,
    variables: Optional[Union[str, int, List[Union[str, int]]]] = None,
    include_flags: bool = True,
    groups: int = 1,
    mf_dataset_freq: Optional[str] = None,
    temp_directory: Optional[Union[str, Path]] = None,
    n_workers: int = 1,
) -> None:
    """

    Parameters
    ----------
    source_files: Union[str, Path]
    output_folder: Union[str, Path]
    variables: Optional[Union[str, int, List[Union[str, int]]]]
    time_step: {"hourly", "daily"}
    include_flags: bool
    groups: int
      The number of file groupings used for converting to multi-file Datasets.
    mf_dataset_freq: Optional[str]
      Resampling frequency for creating output multi-file Datasets. E.g. 'YS': 1 year per file, '5YS': 5 years per file.
    temp_directory: Optional[Union[str, Path]]
      Use another temporary directory location in case default location is not spacious enough.
    n_workers: int

    Returns
    -------
    None
    """
    func_time = time.time()

    if isinstance(source_files, str):
        source_files = Path(source_files)

    if time_step.lower() in ["h", "hour", "hourly"]:
        mode = "hourly"
    elif time_step.lower() in ["d", "day", "daily"]:
        mode = "daily"
    else:
        raise ValueError("Time step must be `h` / `hourly` or `d` / `daily`.")

    if isinstance(variables, list):
        pass
    elif isinstance(variables, (str, int)):
        variables = [variables]
    # TODO: have the variable gathered from a JSON file
    elif variables is None:
        if mode == "hourly":
            variables = [
                89,
                94,
                123,
            ]
            variables.extend(range(76, 81))
            variables.extend(range(262, 281))
        elif mode == "daily":
            variables = [1, 2, 3]
            variables.extend(range(10, 26))
    else:
        raise NotImplementedError()

    for variable_code in variables:
        info = cf_station_metadata(variable_code)
        variable_name = info["nc_name"]
        logging.info(f"Merging `{variable_name}` using `{time_step}` time step.")

        # Only perform aggregation on available data with corresponding metadata
        logging.info("Performing glob and sort.")
        nc_list = list(source_files.joinpath(variable_name).rglob("*.nc"))

        if nc_list != list():
            nc_lists = np.array_split(nc_list, groups)

            with tempfile.TemporaryDirectory(
                prefix="eccc", dir=temp_directory
            ) as temp_dir:

                combinations = sorted(
                    (ii, nc, temp_dir, groups) for ii, nc in enumerate(nc_lists)
                )

                with mp.Pool(processes=n_workers) as pool:
                    pool.starmap(_tmp_zarr, combinations)
                    pool.close()
                    pool.join()

                # TODO memory use seems ok here .. could try using Pool() to increase performance

                zarrs_found = [f for f in Path(temp_dir).glob("*.zarr")]
                logging.info(
                    f"Found {len(zarrs_found)} intermediary aggregation files."
                )

                ds = xr.open_mfdataset(
                    zarrs_found,
                    engine="zarr",
                    combine="nested",
                    concat_dim={"station"},
                )

                if ds:
                    station_file_codes = [x.name.split("_")[0] for x in nc_list]
                    if not include_flags:
                        drop_vars = [vv for vv in ds.data_vars if "flag" in vv]
                        ds = ds.drop_vars(drop_vars)
                    ds = ds.sortby(ds.station_id)

            # Rearrange column order to have lon, lat, elev first
            # # FIXME: This doesn't work as intended - Assign coordinates instead
            # cols = meta.columns.tolist()
            # cols1 = [
            #     "latitude",
            #     "longitude",
            #     "elevation",
            # ]
            # for rr in cols1:
            #     cols.remove(rr)
            # cols1.extend(cols)
            # meta = meta[cols1]
            # meta.index.rename("station", inplace=True)
            # meta = meta.to_xarray()
            # meta.sortby(meta["climate_identifier"])
            # meta = meta.assign({"station": ds.station.values})

            # np.testing.assert_array_equal(
            #     sorted(meta["climate_identifier"].values), sorted(ds.station_id.values)
            # )
            # for vv in meta.data_vars:
            #     ds = ds.assign_coords({vv: meta[vv]})
            # ds = xr.merge([ds, meta])
            # ds.attrs = attrs1

            valid_stations = list(sorted(ds.station_id.values))
            valid_stations_count = len(valid_stations)

            logging.info(f"Processing stations for variable `{variable_name}`.")

            if len(station_file_codes) == 0:
                logging.error(
                    f"No stations were found containing variable filename `{variable_name}`. Exiting."
                )
                return

            logging.info(
                f"Files exist for {len(station_file_codes)} ECCC stations. "
                f"Metadata found for {valid_stations_count} stations. "
            )

            # FIXME: Is this still needed?
            # logging.info("Preparing the NetCDF time period.")
            # Create the time period timestamps
            # year_start = ds.time.dt.year.min().values
            # year_end = ds.time.dt.year.max().values

            # Calculate the time index dimensions of the output NetCDF
            # time_index = pd.date_range(
            #     start=f"{year_start}-01-01",
            #     end=f"{year_end + 1}-01-01",
            #     freq=mode[0].capitalize(),
            # )[:-1]
            # logging.info(
            #     f"Number of ECCC stations: {valid_stations_count}, time steps: {time_index.size}."
            # )

            output_folder.mkdir(parents=True, exist_ok=True)
            file_out = Path(output_folder).joinpath(f"{variable_name}_eccc_{mode}")

            ds = ds.assign_coords(station=range(0, len(ds.station)))
            if mf_dataset_freq is not None:
                # output mf_dataset using resampling frequency
                _, datasets = zip(*ds.resample(time=mf_dataset_freq))
            else:
                datasets = [ds]

            paths = [
                f"{file_out}_{data.time.dt.year.min().values}-{data.time.dt.year.max().values}.nc"
                for data in datasets
            ]

            # FIXME: chunks need to be dealt with
            # chunks = [1, len(ds.time)]
            comp = dict(zlib=True, complevel=5)  # , chunk sizes=chunks)

            with ProgressBar():
                for dataset, path in zip(datasets, paths):
                    encoding = {var: comp for var in ds.data_vars}
                    dataset.to_netcdf(
                        path,
                        engine="netcdf4",
                        format="NETCDF4_CLASSIC",
                        encoding=encoding,
                    )
                    dataset.close()
                    del dataset
            ds.close()

        else:
            logging.info(f"No files found for variable: `{variable_name}`.")

    runtime = f"Process completed in {time.time() - func_time:.2f} seconds"
    logging.warning(runtime)


def _tmp_zarr(
    iterable: int,
    nc: Union[str, Path],
    tempdir: Union[str, Path],
    group: Optional[int] = None,
) -> None:

    logging.info(
        f"Processing batch of files {iterable + 1}"
        f"{' of ' + str(group) if group is not None else ''}."
    )
    station_file_codes = [x.name.split("_")[0] for x in nc]

    ds = xr.open_mfdataset(nc, combine="nested", concat_dim={"station"})
    ds = ds.assign_coords(
        station_id=xr.DataArray(station_file_codes, dims="station").astype(str)
    )
    if "flag" in ds.data_vars:
        ds1 = ds.drop_vars("flag").copy(deep=True)
        ds1["flag"] = ds.flag.astype(str)
        ds = ds1

    with ProgressBar():
        ds.load().to_zarr(
            Path(tempdir).joinpath(f"{str(iterable).zfill(3)}.zarr"),
        )
    del ds


def _combine_years(
    station_folder: str,
    varia: str,
    out_folder: Union[str, Path],
    meta_file: Union[str, Path],
    rejected: List[str],
    _verbose: bool = False,
) -> None:

    nc_files = sorted(list(Path(station_folder).glob("*.nc")))
    logging.info(
        f"Found {len(nc_files)} files for station code {Path(station_folder).name}."
    )

    # Remove range files if years are all present, otherwise default to range_file.
    years_found = dict()
    range_files_found = dict()
    years_parsed = True
    for f in nc_files:
        groups = re.findall(r"_\d{4}", f.stem)
        if len(groups) == 1:
            year = int(groups[0].strip("_"))
            years_found[year] = f
        elif len(groups) == 2:
            year_start, year_end = int(groups[0].strip("_")), int(groups[1].strip("_"))
            range_files_found[f] = set(range(year_start, year_end))
        else:
            logging.warning(
                "Years unable to be effectively parsed from series. Continuing with xarray solver..."
            )
            years_parsed = False
            break
    if years_parsed:
        if len(range_files_found) > 0:
            logging.warning(
                f"Overlapping single-year and multi-year files found for {station_folder}. Removing overlaps."
            )
            for ranged_file, years in range_files_found.items():
                if years.issubset(years_found.values()):
                    nc_files.remove(ranged_file)
                else:
                    for y in years:
                        try:
                            nc_files.remove(years_found[y])
                        except (KeyError, ValueError):
                            continue

        year_range = min(years_found.keys()), max(years_found.keys())
        logging.info(
            "Year(s) covered: "
            f"{year_range[0]}{'-' + str(year_range[1]) if year_range[0] != year_range[1] else ''}. "
        )

    if _verbose:
        logging.info(f"Opening: {', '.join([p.name for p in nc_files])}")
    ds = xr.open_mfdataset(nc_files, combine="nested", concat_dim={"time"})
    outfile = out_folder.joinpath(
        f'{nc_files[0].name.split(f"_{varia}_")[0]}_{varia}_'
        f"{ds.time.dt.year.min().values}-{ds.time.dt.year.max().values}.nc"
    )

    df_inv = xr.open_dataset(meta_file)

    station_id = ds.attrs["member"]
    meta = df_inv.isel(index=df_inv.CLIMATE_IDENTIFIER == station_id)
    meta = meta.rename({"index": "station", "CLIMATE_IDENTIFIER": "station_id"})
    try:
        meta = meta.assign_coords(station=[0])
    except ValueError:
        rejected.append(Path(station_folder).name)
        logging.error(
            f"Something went wrong at the assign_coords step for station {station_folder}. Continuing..."
        )
        return
    if len(meta.indexes) > 1:
        raise ValueError("Found more than 1 station.")
    elif len(meta.indexes) == 0:
        rejected.append(Path(station_folder).name)
        logging.warning(
            f"No metadata found for station {station_folder}. Continuing..."
        )
        return

    keep_coords = [
        "time",
        "station",
        "station_id",
        "latitude",
        "longitude",
        "elevation",
    ]
    for vv in meta.data_vars:
        if vv.lower() not in keep_coords:
            continue
        ds = ds.assign_coords({vv.lower(): meta[vv]})

    for vv in ds.data_vars:
        if ds[vv].dtype == "O":
            ds[vv] = ds[vv].astype(str)

    if not outfile.exists():
        logging.info(f"Merging to {outfile.name}")
        comp = dict(zlib=True, complevel=5)
        encoding = {data_var: comp for data_var in ds.data_vars}
        encoding["time"] = dict(dtype="single")
        with ProgressBar():
            ds.to_netcdf(
                outfile,
                engine="netcdf4",
                format="NETCDF4",
                encoding=encoding,
            )
    else:
        logging.info(f"Files exist for {outfile.name}. Continuing...")


def merge_converted_variables(
    source: Union[str, Path],
    destination: Union[str, Path],
    variables: Optional[Union[str, int, List[Union[str, int]]]] = None,
    station_metadata: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    n_workers: int = 1,
) -> None:
    """

    Parameters
    ----------
    source: Union[str, Path]
    destination: Union[str, Path]
    variables: Optional[Union[str, int, List[Union[str, int]]]]
    station_metadata: Optional[Union[str, Path]]
    overwrite: bool
    n_workers: int

    Returns
    -------
    None

    """
    meta = load_station_metadata(station_metadata)
    metadata_file = Path(tempfile.NamedTemporaryFile(suffix=".nc", delete=False).name)
    meta.to_netcdf(metadata_file)

    if isinstance(source, str):
        source = Path(source)
    if isinstance(destination, str):
        destination = Path(destination)

    selected_variables = list()
    if variables is not None:
        if not isinstance(variables, list):
            variables = [variables]
        for var in variables:
            selected_variables.append(cf_station_metadata(var))

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

        output_rep = destination.joinpath(variable)
        Path(output_rep).mkdir(parents=True, exist_ok=True)

        if (
            len(list(output_rep.iterdir())) >= (len(meta.CLIMATE_IDENTIFIER) * 0.75)
        ) and not overwrite:
            logging.warning(
                f"Variable {variable} appears to have already been converted. Will be skipped. "
                f"To force conversion of this variable, set `overwrite=True`."
            )
            continue

        manager = mp.Manager()
        rejected_stations = manager.list()

        combine_func = functools.partial(
            _combine_years,
            varia=variable,
            out_folder=output_rep,
            meta_file=metadata_file,
            rejected=rejected_stations,
        )

        with mp.Pool(processes=n_workers) as pool:
            pool.map(combine_func, station_dirs)
            pool.close()
            pool.join()

        if rejected_stations:
            logging.warning(
                f"Rejected station codes are the following: {', '.join(rejected_stations)}."
            )
