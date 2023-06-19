######################################################################
# G. Rondeau-Genesse, Ouranos, 2019-09-27
#
# Description
#
# find_and_extract_dly finds all the CSV files of a ECCC daily weather station,
# then appends the data within a pandas Dataframe
#
# dly_to_netcdf takes that Dataframe and exports it to a netCDF. When possible,
# the variables are converted to be compatible with CF-Convention. For example,
# "Max Temp (°C)" is renamed "tasmax" and converted to °K.
#
#####################################################################
from __future__ import annotations

import json
import logging
from collections import defaultdict
from collections.abc import Generator
from logging import config
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from miranda.scripting import LOGGING_CONFIG

config.dictConfig(LOGGING_CONFIG)
__all__ = ["extract_daily_summaries", "daily_summaries_to_netcdf"]

eccc_metadata = json.load(
    open(Path(__file__).parent / "eccc_obs_summary_cf_attrs.json")
)["variable_entry"]


# Searches a location for the station data, then calls the needed scripts to read and assembles the data using pandas
def extract_daily_summaries(
    path_station: Path | str, rm_flags: bool = False, file_suffix: str = ".csv"
) -> dict:
    """Extract daily climate summaries from ECCC CSV files.

    Parameters
    ----------
    path_station : str or Path
        PathLike or str to the station's folder containing the csv files.
    rm_flags : bool
        Removes the 'Flag' and 'Quality' columns of the ECCC files.
    file_suffix : str
        File suffixes used by the tabular data. Default: ".csv".

    Returns
    -------
    dict
        dict containing the station metadata, as well as the data stored within a pandas Dataframe.
    """
    # Find the CSV files
    if "*" not in file_suffix:
        file_suffix = f"*{file_suffix}"
    station_files = Path(path_station).rglob(file_suffix)

    # extract the .csv data
    stations = _read_multiple_daily_summaries(station_files, rm_flags=rm_flags)

    return stations


#
def daily_summaries_to_netcdf(station: dict, path_output: Path | str) -> None:
    """Convert daily climate summaries to NetCDF files.

    Uses xarray to transform the 'station' from find_and_extract_dly into a CF-Convention netCDF file

    Parameters
    ----------
    station : dict
        dict created by using find_and_extract_dly
    path_output: str or Path
        Output path.

    Returns
    -------
    None
    """
    # first, transform the Date/Time to a 'days since' format
    time = station["data"]["Date/Time"] - np.array(
        "1950-01-01T00:00", dtype="datetime64"
    )
    time = time.astype("timedelta64[s]").astype(float) / 86400

    # we use expand_dims twice to 'add' longitude and latitude dimensions to the station data
    logging.info(
        "Reading data for station {} (ID: {}) now.".format(
            station["name"], station["ID"]
        )
    )

    ds = None

    variables = eccc_metadata["variable_entry"]
    for var in variables.keys():
        original_field = variables[var]["original_field"]
        add_offset = variables[var]["add_offset"]
        scale_factor = variables[var]["scale_factor"]

        da = xr.DataArray(
            np.expand_dims(
                np.expand_dims(
                    station["data"][original_field] * scale_factor + add_offset, axis=1
                ),
                axis=2,
            ),
            [
                ("time", time),
                ("lat", [station["latitude"]]),
                ("lon", [station["longitude"]]),
            ],
        )

        da.name = var
        for field in [
            "standard_name",
            "long_name",
            "units",
            "grid_mapping",
            "comments",
            "frequency",
        ]:
            da.attrs[field] = variables[var][field]

        # for the first variable, we simply create a dataset from it
        if ds is None:
            ds = da.to_dataset()
        else:
            ds[var] = da

    # add attributes to lon, lat, time, elevation, and the grid
    # TODO: There is probably a better CF Convention for point-based data
    da = xr.DataArray(np.full(len(time), np.nan), [("time", time)])
    da.name = "regular_lon_lat"
    da.attrs["grid_mapping_name"] = "lonlat"
    ds["regular_lon_lat"] = da

    da = xr.DataArray(
        np.expand_dims(np.expand_dims(station["elevation"], axis=1), axis=2),
        [("lat", [station["latitude"]]), ("lon", [station["longitude"]])],
    )
    da.name = "elevation"
    da.attrs["standard_name"] = "elevation"
    da.attrs["long_name"] = "elevation"
    da.attrs["units"] = "m"
    da.attrs["axis"] = "Z"
    ds["elevation"] = da
    ds = ds.set_coords("elevation")

    ds.lon.attrs["standard_name"] = "longitude"
    ds.lon.attrs["long_name"] = "longitude"
    ds.lon.attrs["units"] = "degrees_east"
    ds.lon.attrs["axis"] = "X"

    ds.lat.attrs["standard_name"] = "latitude"
    ds.lat.attrs["long_name"] = "latitude"
    ds.lat.attrs["units"] = "degrees_north"
    ds.lat.attrs["axis"] = "Y"

    ds.time.attrs["standard_name"] = "time"
    ds.time.attrs["long_name"] = "time"
    ds.time.attrs["units"] = "days since 1950-01-01 00:00:00"
    ds.time.attrs["axis"] = "T"
    ds.time.attrs["calendar"] = "gregorian"

    # add global attributes
    ds.attrs["Station Name"] = station["name"]
    ds.attrs["Province"] = station["province"]
    ds.attrs["Climate Identifier"] = station["ID"]
    ds.attrs["WMO Identifier"] = station["WMO_ID"]
    ds.attrs["TC Identifier"] = station["TC_ID"]
    ds.attrs["Institution"] = "Environment and Climate Change Canada"

    # save the data
    output_file = Path(path_output).joinpath("{}.nc".format(ds.attrs["Station Name"]))
    ds.to_netcdf(output_file)


##########################################
# BELOW THIS POINT ARE UTILITY SCRIPTS
##########################################


# This
def _read_multiple_daily_summaries(
    files: list[str | Path] | Generator[Path, None, None],
    rm_flags: bool = False,
) -> dict:
    """

    Notes
    -----
    This calls `_read_single_eccc_dly` and appends the data in a single Dict.

    Parameters
    ----------
    files : list of str or Path, or Generator[Path]
        A list of all the files to append.
    rm_flags : bool
        Removes all the 'Flag' and 'Quality' columns of the ECCC files. Default: False.

    Returns
    -------
    dict
    """
    # Extract the data for each files
    all_stations = dict()
    station_data = list()

    file_list = [Path(f) for f in files]
    file_list.sort()

    station_codes = defaultdict(list)
    for f in file_list:
        code = Path(f).name.split("_")[4]
        station_codes[code].append(f)

    for station_code, summary_files in station_codes.items():
        for summary in summary_files:
            station = pd.read_csv(summary)
            station_data.append(station)

        station_summary_full = pd.DataFrame(
            station_data
        )  # FIXME: Find the way to combine list of dataframes into one

        # change the Date/Time column to a datetime64 type
        station_summary_full["Date/Time"] = pd.to_datetime(
            station_summary_full["Date/Time"]
        )

        # if wanted, remove the quality and flag columns
        if rm_flags:
            index_quality = [
                i
                for i, s in enumerate(station_summary_full.columns.values)
                if "Quality" in s
            ]
            station_summary_full = station_summary_full.drop(
                station_summary_full.columns.values[index_quality], axis="columns"
            )
            index_flag = [
                i
                for i, s in enumerate(station_summary_full.columns.values)
                if "Flag" in s
            ]
            station_summary_full = station_summary_full.drop(
                station_summary_full.columns.values[index_flag], axis="columns"
            )

        # combine everything in a single Dict
        all_stations[station_code] = station_summary_full

    return all_stations


def _read_single_daily_summaries(file: str | Path) -> tuple[dict, pd.DataFrame]:
    """Read station summary information from CSV header.

    Notes
    -----
    Climate Services Canada has changed the way they store metadata and no longer store this infor in the CSV heading.

    Parameters
    ----------
    file : str or Path

    Returns
    -------
    tuple[dict, pd.DataFrame]
    """
    # Read the whole file
    with open(file, encoding="utf-8-sig") as fi:
        lines = fi.readlines()

    # Find each element in the header
    search_header = [0] * 9
    search_header[0] = [i for i, s in enumerate(lines) if "Station Name" in s][0]
    search_header[1] = [i for i, s in enumerate(lines) if "Province" in s][0]
    search_header[2] = [i for i, s in enumerate(lines) if "Latitude" in s][0]
    search_header[3] = [i for i, s in enumerate(lines) if "Longitude" in s][0]
    search_header[4] = [i for i, s in enumerate(lines) if "Elevation" in s][0]
    search_header[5] = [i for i, s in enumerate(lines) if "Climate Identifier" in s][0]
    search_header[6] = [i for i, s in enumerate(lines) if "WMO Identifier" in s][0]
    search_header[7] = [i for i, s in enumerate(lines) if "TC Identifier" in s][0]
    search_header[8] = [i for i, s in enumerate(lines) if "Date/Time" in s][0]
    # This is where the data actually starts

    # Does a bunch of stuff, but basically finds the right line, then cleans up the string
    station_meta = {
        "name": lines[search_header[0]]
        .split(",")[1]
        .replace('"', "")
        .replace("\n", ""),
        "province": lines[search_header[1]]
        .split(",")[1]
        .replace('"', "")
        .replace("\n", ""),
        "latitude": float(
            lines[search_header[2]].split(",")[1].replace('"', "").replace("\n", "")
        ),
        "longitude": float(
            lines[search_header[3]].split(",")[1].replace('"', "").replace("\n", "")
        ),
        "elevation": float(
            lines[search_header[4]].split(",")[1].replace('"', "").replace("\n", "")
        ),
        "ID": lines[search_header[5]].split(",")[1].replace('"', "").replace("\n", ""),
        "WMO_ID": lines[search_header[6]]
        .split(",")[1]
        .replace('"', "")
        .replace("\n", ""),
        "TC_ID": lines[search_header[7]]
        .split(",")[1]
        .replace('"', "")
        .replace("\n", ""),
    }

    data = pd.read_csv(file, header=search_header[8] - 2)
    # Makes sure that the data starts on Jan 1st
    if data.values[0, 2] != 1 | data.values[0, 3] != 1:
        logging.warning(
            f"Data for file {file.name} is not starting on January 1st. Make sure this is what you want!"
        )

    return station_meta, data
