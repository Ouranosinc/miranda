"""Specialized conversion tools for Environment and Climate Change Canada / Meteorological Service of Canada data."""

from __future__ import annotations
import functools
import logging
import multiprocessing as mp
import os
import re
import tempfile
import time
from datetime import datetime as dt
from pathlib import Path
from typing import Any

import dask.dataframe as dd
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar

from miranda.preprocess._metadata import eccc_variable_metadata, obs_column_definitions
from miranda.treatments import find_project_variable_codes, load_json_data_mappings
from miranda.vocabularies.eccc import obs_vocabularies


__all__ = [
    "convert_station",
    "merge_converted_variables",
    "merge_stations",
]
TABLE_DATE = dt.now().strftime("%d %B %Y")


def _remove_duplicates(ds):
    if any(ds.get_index("time").duplicated()):
        msg = f"Found {ds.get_index('time').duplicated().sum()} duplicated time coordinates for station {ds.station_id.values}. Assuming first value."
        logging.info(msg)
    return ds.sel(time=~ds.get_index("time").duplicated())


def convert_observation(
    data_source: str | Path | list[str | Path],
    output_dir: str | Path,
    variable: str,
    *,
    generation: int | None = None,
    merge: bool = False,
    overwrite: bool = False,
):
    """Convert a single station's data from the fixed-width format to a netCDF file."""
    output_dir = Path(output_dir).resolve().joinpath(variable)
    output_dir.mkdir(parents=True, exist_ok=True)

    code = find_project_variable_codes(variable, "eccc-obs")
    var_meta, global_attrs = eccc_variable_metadata(code, "eccc-obs", generation)
    (
        column_names,
        column_spaces,
        column_dtypes,
        header_row,
    ) = obs_column_definitions(code)

    archives = list()
    if isinstance(data_source, list) or Path(data_source).is_file():
        archives.append(data_source)
    else:
        tables = [str(repository.keys()) for repository in obs_vocabularies if code in repository.values()]
        msg = f"Collecting files for variable '{variable}'. Filename patterns containing variable code '{code}: {', '.join(tables)}'."
        logging.info(msg)
        for table in tables:
            archives.extend([f for f in Path(data_source).rglob(f"{table}*.gz")])

    # Create the output directory
    output_variable_dir = Path(output_dir).joinpath(variable)
    output_variable_dir.mkdir(parents=True, exist_ok=True)

    # Loop on the files
    errored_files = []
    for _file in archives:
        # FIXME: convert the file using the appropriate function
        pass

    if errored_files:
        msg = "Some files failed to be properly parsed:\n", ", ".join(errored_files)
        logging.warning(msg)


def convert_station(
    data: str | os.PathLike,
    variable: str,
    mode: str,
    # missing_flags: set[str],
    # missing_values: set[str],
    using_dask_array: bool = False,
    *,
    client: Any,
    **kwargs,
):
    """Convert a single station's data from the fixed-width format to a netCDF file."""
    data = Path(data)
    variable_code = find_project_variable_codes(variable, "eccc-obs")
    column_names, column_widths, column_dtypes, header = obs_column_definitions(mode)

    # if not missing_values:
    #     missing_values = {-9999, "#####"}

    if using_dask_array:
        pandas_reader = dd
        # set the block size to 200 MB
        chunks = dict(blocksize=200 * 2**20)
    else:
        pandas_reader = pd
        chunks = dict()
        using_dask_array = False

    # Create a dataframe from the files
    try:
        df = pandas_reader.read_fwf(
            data,
            widths=column_widths,
            names=column_names,
            dtype={name: data_type for name, data_type in zip(column_names, column_dtypes, strict=False)},
            assume_missing=True,
            **chunks,
        )
        if using_dask_array:
            df = client.persist(df)

    except FileNotFoundError as err:
        msg = f"File {data} was not found."
        logging.error(msg)
        raise FileNotFoundError(msg) from err

    except UnicodeDecodeError as e:
        msg = f"File {data.name} was unable to be read. This is probably an issue with the file: {e}"
        logging.error(msg)
        raise

    # Loop through the station codes
    station_codes = df["code"].unique()
    for code in station_codes:
        df_code = df[df["code"] == code]

        # Abort if the variable is not found
        if using_dask_array:
            has_variable_codes = ((df_code["code_var"] == variable_code).compute()).any()
        else:
            has_variable_codes = (df_code["code_var"] == variable_code).any()
        if not has_variable_codes:
            msg = f"Variable `{variable}` not found for station code: {code} in file {data}. Continuing..."
            logging.info(msg)
            continue

        # # Perform the data treatment
        # logging.info(f"Converting `{variable}` for station code: {code}")
        #
        # # Dump the data into a DataFrame
        # df_var = df_code[df_code["code_var"] == variable_code].copy()
        #
        # # Mask the data according to the missing values flag
        # df_var = df_var.replace(missing_values, np.nan)
        #
        # # Decode the values and flags
        # dfd = df_var.loc[:, [f"D{i:0n}" for i in range(1, num_observations + 1)]]
        # dff = df_var.loc[:, [f"F{i:0n}" for i in range(1, num_observations + 1)]]
        #
        # # Remove the "NaN" flag
        # dff = dff.fillna("")
        #
        # # Use the flag to mask the values
        # try:
        #     val = np.asarray(dfd.values, float)
        # except ValueError as e:
        #     logging.error(f"{e} raised from {dfd}, continuing...")
        #     continue
        # try:
        #     flag = np.asarray(dff.values, str)
        # except ValueError as e:
        #     logging.error(f"{e} raised from {dff}, continuing...")
        #     continue
        # mask = np.isin(flag, missing_flags)
        # val[mask] = np.nan
        #
        # # Treat according to units conversions
        # val = val * scale_factor + add_offset

        # Create the DataArray
        # date_summations = dict(time=list())
        # if mode == "hourly":
        #     for index, row in df_var.iterrows():
        #         period = pd.Period(
        #             year=row.year, month=row.month, day=row.day, freq="D"
        #         )
        #         dates = pd.Series(
        #             pd.date_range(
        #                 start=period.start_time,
        #                 end=period.end_time,
        #                 freq="H",
        #             )
        #         )
        #         date_summations["time"].extend(dates)
        #     written_values = val.flatten()
        #     written_flags = flag.flatten()
    #     elif mode == "daily":
    #         value_days = list()
    #         flag_days = list()
    #         for i, (index, row) in enumerate(df_var.iterrows()):
    #             period = pd.Period(year=row.year, month=row.month, freq="M")
    #             dates = pd.Series(
    #                 pd.date_range(
    #                     start=period.start_time,
    #                     end=period.end_time,
    #                     freq="D",
    #                 )
    #             )
    #             date_summations["time"].extend(dates)
    #
    #             value_days.extend(
    #                 val[i][range(monthrange(int(row.year), int(row.month))[1])]
    #             )
    #             flag_days.extend(
    #                 flag[i][range(monthrange(int(row.year), int(row.month))[1])]
    #             )
    #         written_values = value_days
    #         written_flags = flag_days
    #
    #     ds = xr.Dataset()
    #     da_val = xr.DataArray(written_values, coords=date_summations, dims=["time"])
    #
    #     if raw_units != units:
    #         da_val.attrs["units"] = raw_units
    #         da_val = convert_units_to(da_val, units)
    #     else:
    #         da_val.attrs["units"] = units
    #
    #     da_val = da_val.rename(nc_name)
    #     variable_attributes = dict(
    #         variable_code=variable_code,
    #         standard_name=standard_name,
    #         long_name=long_name,
    #     )
    #     if "original_units" in kwargs:
    #         variable_attributes["original_units"] = kwargs["original_units"]
    #     da_val.attrs.update(variable_attributes)
    #
    #     da_flag = xr.DataArray(written_flags, coords=date_summations, dims=["time"])
    #     da_flag = da_flag.rename("flag")
    #     flag_attributes = dict(
    #         long_name="data flag",
    #         note="See ECCC technical documentation for details",
    #     )
    #     da_flag.attrs.update(flag_attributes)
    #
    #     ds[nc_name] = da_val
    #     ds["flag"] = da_flag
    #
    #     # save the file in NetCDF format
    #     start_year = ds.time.dt.year.values[0]
    #     end_year = ds.time.dt.year.values[-1]
    #
    #     station_folder = output_path.joinpath(str(code))
    #     station_folder.mkdir(parents=True, exist_ok=True)
    #
    #     f_nc = (
    #         f"{code}_{variable_code}_{nc_name}_"
    #         f"{start_year if start_year == end_year else '_'.join([str(start_year), str(end_year)])}.nc"
    #     )
    #
    #     if station_folder.joinpath(f_nc).exists():
    #         logging.warning(f"File `{f_nc}` already exists. Continuing...")
    #
    #     history = (
    #         f"{dt.now().strftime('%Y-%m-%d %X')} converted from flat station file "
    #         f"(`{file.name}`) to n-dimensional array."
    #     )
    #
    #     # TODO: This info should eventually be sourced from a JSON definition
    #     global_attrs = dict(
    #         Conventions="CF-1.8",
    #         comment="Acquired on demand from data specialists at "
    #         "ECCC Climate Services / Services Climatiques.",
    #         contact="John Richard",
    #         contact_email="climatcentre-climatecentral@ec.gc.ca",
    #         domain="CAN",
    #     )
    #     if mode == "hourly":
    #         global_attrs.update(dict(frequency="1hr"))
    #     elif mode == "daily":
    #         global_attrs.update(dict(frequency="day"))
    #     global_attrs.update(
    #         dict(
    #             history=history,
    #             internal_comment=f"Converted by {os.environ.get('USER', os.environ.get('USERNAME'))}.",
    #             institution="ECCC",
    #             license="https://climate.weather.gc.ca/prods_servs/attachment1_e.html",
    #             member=code,
    #             processing_level="raw",
    #             redistribution="Redistribution permitted.",
    #             references="https://climate.weather.gc.ca/doc/Technical_Documentation.pdf",
    #             source="historical-station-records",
    #             table_date=TABLE_DATE,
    #             title="Environment and Climate Change Canada (ECCC) weather station observations",
    #             type="station-obs",
    #             usage="The original data is owned by the Government of Canada (Environment and Climate "
    #             "Change Canada), and falls under the licence agreement for use of Environment and "
    #             "Climate Change Canada data",
    #             variable=str(nc_name),
    #             version=f"v{dt.now().strftime('%Y.%m.%V')}",  # Year.Month.Week
    #         )
    #     )
    #     ds.attrs.update(global_attrs)
    #
    #     logging.info(f"Exporting to: {station_folder.joinpath(f_nc)}")
    #     ds.to_netcdf(station_folder.joinpath(f_nc))
    #     del ds
    #     del val
    #     del mask
    #     del flag
    #     del da_val
    #     del da_flag
    #     del dfd
    #     del dff
    #     del written_values
    #     del written_flags
    #     del date_summations
    #
    # del df


def merge_stations(
    source_files: str | os.PathLike | None = None,
    output_folder: str | os.PathLike | None = None,
    *,
    time_step: str,
    variables: str | int | list[str | int] | None = None,
    include_flags: bool = True,
    groupings: int | None = None,
    mf_dataset_freq: str | None = None,
    temp_directory: str | os.PathLike | None = None,
    n_workers: int = 1,
) -> None:
    """
    Merge stations.

    Parameters
    ----------
    source_files : str or Path
        Source files to be aggregated.
    output_folder : str or Path
        Output folder for the aggregated files.
    variables : str or int or list of str or int, optional
        The variable codes to be aggregated.
    time_step : {"hourly", "daily"}
        The time step to be used for aggregation.
    include_flags : bool
        Include flags in the output files.
    groupings : int
        The number of files in each group used for converting to multi-file Datasets.
    mf_dataset_freq : str, optional
        Resampling frequency for creating output multi-file Datasets. E.g. 'YS': 1 year per file, '5YS': 5 years per file.
    temp_directory : str or Path, optional
        Use another temporary directory location in case default location is not spacious enough.
    n_workers : int
        The number of workers to use.

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
        info = load_json_data_mappings("eccc-obs")["variables"][variable_code]
        variable_name = info["cf_variable_name"]
        msg = f"Merging `{variable_name}` using `{time_step}` time step."
        logging.info(msg)

        # Only perform aggregation on available data with corresponding metadata
        logging.info("Performing glob and sort.")
        nc_list = [str(nc) for nc in source_files.joinpath(variable_name).rglob("*.nc")]

        if not groupings:
            groupings = max(n_workers**2, 4)

        if nc_list:
            with tempfile.TemporaryDirectory(prefix="eccc", dir=temp_directory) as temp_dir:
                combinations = sorted((ii, nc, temp_dir, len(nc_list)) for ii, nc in enumerate(nc_list))

                with mp.Pool(processes=n_workers) as pool:
                    pool.starmap(_tmp_zarr, combinations)
                    pool.close()
                    pool.join()

                zarrs_found = [f for f in Path(temp_dir).glob("*.zarr")]
                msg = f"Found {len(zarrs_found)} intermediary aggregation files."
                logging.info(msg)

                ds = xr.open_mfdataset(
                    zarrs_found,
                    engine="zarr",
                    combine="nested",
                    concat_dim={"station"},
                )

                if ds:
                    station_file_codes = [Path(x).name.split("_")[0] for x in nc_list]
                    if not include_flags:
                        drop_vars = [vv for vv in ds.data_vars if "flag" in vv]
                        ds = ds.drop_vars(drop_vars)
                    ds = ds.sortby(ds.station_id, "time")

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

                # export done within tmddir context otherwise data is erased before final export!!
                valid_stations = list(sorted(ds.station_id.values))
                valid_stations_count = len(valid_stations)

                msg = f"Processing stations for variable `{variable_name}`."
                logging.info(msg)

                if len(station_file_codes) == 0:
                    msg = f"No stations were found containing variable filename `{variable_name}`. Exiting."
                    logging.error(msg)
                    return

                msg = f"Files exist for {len(station_file_codes)} ECCC stations. Metadata found for {valid_stations_count} stations. "
                logging.info(msg)

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

                Path(output_folder).mkdir(parents=True, exist_ok=True)
                file_out = Path(output_folder).joinpath(f"{variable_name}_eccc_{mode}")

                ds = ds.assign_coords(station=range(0, len(ds.station))).sortby("time")
                if mf_dataset_freq is not None:
                    # output mf_dataset using resampling frequency
                    _, datasets = zip(*ds.resample(time=mf_dataset_freq), strict=False)
                else:
                    datasets = [ds]

                paths = [f"{file_out}_{data.time.dt.year.min().values}-{data.time.dt.year.max().values}.nc" for data in datasets]

                # FIXME: chunks need to be dealt with
                # chunks = [1, len(ds.time)]
                # comp = dict(zlib=True, complevel=5)  # , chunk sizes=chunks)

                with ProgressBar():
                    # FIXME: looping seems to cause increasing memory over time use a pool of one or 2??
                    # for dataset, path in zip(datasets, paths):
                    #     _export_agg_nc(dataset,path)
                    combs = zip(datasets, paths, strict=False)
                    pool = mp.Pool(2)
                    pool.map(_export_agg_nc, combs)
                    pool.close()
                    pool.join()
                ds.close()
                del ds

        else:
            msg = f"No files found for variable: `{variable_name}`."
            logging.info(msg)

    runtime = f"Process completed in {time.time() - func_time:.2f} seconds."
    logging.warning(runtime)


def _export_agg_nc(args):
    dataset, path = args
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in dataset.data_vars}
    dataset.load().to_netcdf(
        path,
        engine="h5netcdf",
        format="NETCDF4_CLASSIC",
        encoding=encoding,
    )
    dataset.close()
    del dataset


def _tmp_zarr(
    iterable: int,
    nc: list[str | os.PathLike],
    tempdir: str | os.PathLike,
    group: int | None = None,
) -> None:
    msg = f"Processing batch of files {iterable + 1}{' of ' + str(group) if group is not None else ''}."
    logging.info(msg)
    station_file_codes = [Path(x).name.split("_")[0] for x in nc]

    try:
        ds = xr.open_mfdataset(nc, combine="nested", concat_dim="station", preprocess=_remove_duplicates)
    except ValueError as e:
        errored_nc_files = ", ".join([Path(f).name for f in nc])
        msg = f"Issues found with the following files: [{errored_nc_files}]: {e}"
        logging.error(msg)
        return

    ds = ds.assign_coords(station_id=xr.DataArray(station_file_codes, dims="station").astype(str))
    if "flag" in ds.data_vars:
        ds1 = ds.drop_vars("flag").copy(deep=True)
        ds1["flag"] = ds.flag.astype(str)
        ds = ds1

    with ProgressBar():
        ds.load().to_zarr(
            Path(tempdir).joinpath(f"{str(iterable).zfill(4)}.zarr"),
        )
    del ds


def _combine_years(
    station_folder: str,
    varia: str,
    out_folder: str | os.PathLike,
    meta_file: str | os.PathLike,
    rejected: list[str],
    _verbose: bool = False,
) -> None:
    nc_files = sorted(list(Path(station_folder).glob("*.nc")))
    if len(nc_files):
        msg = f"Found {len(nc_files)} files for station code {Path(station_folder).name}."
        logging.info(msg)
    else:
        msg = f"No readings found for station code {Path(station_folder).name}. Continuing..."
        logging.warning(msg)
        return

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
            msg = "Years unable to be effectively parsed from series. Continuing with xarray solver..."
            logging.warning(msg)
            years_parsed = False
            break
    if years_parsed:
        if len(range_files_found) > 0:
            msg = f"Overlapping single-year and multi-year files found for station code {station_folder}. Removing overlaps."
            logging.warning(msg)
            for ranged_file, years in range_files_found.items():
                if years.issubset(years_found.values()):
                    nc_files.remove(ranged_file)
                else:
                    missing_years = []
                    for y in years:
                        try:
                            nc_files.remove(years_found[y])
                        except (KeyError, ValueError):  # noqa: PERF203
                            missing_years.append(str(y))
                            continue
                        if missing_years:
                            msg = f"Missing years {', '.join(missing_years)} from multi-year file {ranged_file}. "
                            logging.warning(msg)

        year_range = min(years_found.keys()), max(years_found.keys())
        msg = f"Year(s) covered: {year_range[0]}{'-' + str(year_range[1]) if year_range[0] != year_range[1] else ''}. "
        logging.info(msg)

    if _verbose:
        msg = f"Opening: {', '.join([p.name for p in nc_files])}"
        logging.info(msg)
    ds = xr.open_mfdataset(nc_files, combine="nested", concat_dim="time")
    outfile = Path(out_folder).joinpath(
        f"{nc_files[0].name.split(f'_{varia}_')[0]}_{varia}_{ds.time.dt.year.min().values}-{ds.time.dt.year.max().values}.nc"
    )

    df_inv = xr.open_dataset(meta_file)

    station_id = ds.attrs["member"]
    meta = df_inv.isel(index=df_inv.CLIMATE_IDENTIFIER == station_id)
    meta = meta.rename({"index": "station", "CLIMATE_IDENTIFIER": "station_id"})
    try:
        meta = meta.assign_coords(station=[0])
    except ValueError:
        rejected.append(Path(station_folder).name)
        msg = f"Something went wrong at the assign_coords step for station {station_folder}. Continuing..."
        logging.error(msg)
        return
    if len(meta.indexes) > 1:
        raise ValueError("Found more than 1 station.")
    elif len(meta.indexes) == 0:
        rejected.append(Path(station_folder).name)
        msg = f"No metadata found for station code {station_folder}. Continuing..."
        logging.warning(msg)
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
        if str(vv).lower() not in keep_coords:
            continue
        ds = ds.assign_coords({str(vv).lower(): meta[vv]})

    for vv in ds.data_vars:
        if ds[vv].dtype == "O":
            ds[vv] = ds[vv].astype(str)

    if not outfile.exists():
        msg = f"Merging to {outfile.name}."
        logging.info(msg)
        comp = dict(zlib=True, complevel=5)
        encoding = {data_var: comp for data_var in ds.data_vars}
        encoding["time"] = {"dtype": "single"}
        with ProgressBar():
            ds.to_netcdf(
                outfile,
                engine="h5netcdf",
                format="NETCDF4_CLASSIC",
                encoding=encoding,
            )
    else:
        msg = f"Files exist for {outfile.name}. Continuing..."
        logging.info(msg)


def merge_converted_variables(
    source_files: str | os.PathLike,
    output_folder: str | os.PathLike,
    variables: str | int | list[str | int] | None = None,
    overwrite: bool = False,
    n_workers: int = 1,
) -> None:
    """
    Merge converted variables into a single file per variable.

    Parameters
    ----------
    source_files : str, Path
    output_folder : str, Path
    variables : str or int or list of str or int, optional
    overwrite : bool
    n_workers : int

    Returns
    -------
    None
    """
    meta = load_json_data_mappings("eccc-obs")
    metadata_file = Path(tempfile.NamedTemporaryFile(suffix=".nc", delete=False).name)
    meta.to_netcdf(metadata_file)

    if isinstance(source_files, str):
        source_files = Path(source_files)
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

    selected_variables = list()
    if variables is not None:
        if not isinstance(variables, list):
            variables = [variables]
        selected_variables.extend(meta[var] for var in variables)

    variables_found = [x.name for x in source_files.iterdir() if x.is_dir()]
    if selected_variables:
        variables_found = [x for x in variables_found if x in [item["nc_name"] for item in selected_variables]]

    for variable in variables_found:
        msg = f"Merging files found for variable: `{variable}`."
        logging.info(msg)
        station_dirs = [x for x in source_files.joinpath(variable).iterdir() if x.is_dir()]
        msg = f"Number of stations found: {len(station_dirs)}."
        logging.info(msg)

        output_rep = output_folder.joinpath(variable)
        Path(output_rep).mkdir(parents=True, exist_ok=True)

        if (len(list(output_rep.iterdir())) >= (len(meta.CLIMATE_IDENTIFIER) * 0.75)) and not overwrite:
            msg = (
                f"Variable {variable} appears to have already been converted. Will be skipped. "
                f"To force conversion of this variable, set `overwrite=True`."
            )
            logging.warning(msg)
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
            msg = f"Rejected station codes are the following: {', '.join(rejected_stations)}."
            logging.warning(msg)
