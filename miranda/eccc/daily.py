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
import logging
from pathlib import Path
from types import GeneratorType
from typing import Generator
from typing import List
from typing import Tuple
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr

__all__ = ["find_and_extract_dly", "dly_to_netcdf"]


# Searches a location for the station data, then calls the needed scripts to read and assembles the data using pandas
def find_and_extract_dly(
    path_station: Union[Path, str], rm_flags: bool = False, file_suffix: str = ".csv"
) -> dict:
    """

    Parameters
    ----------
    path_station : Union[Path, str]
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
        file_suffix = "*{}".format(file_suffix)
    station_files = Path(path_station).rglob(file_suffix)

    # extract the .csv data
    station = _read_multiple_eccc_dly(station_files, rm_flags=rm_flags)

    return station


# Uses xarray to transform the 'station' from find_and_extract_dly into a CF-Convention netCDF file
def dly_to_netcdf(station: dict, path_output: Union[Path, str]) -> None:
    """

    Parameters
    ----------
    station : dict
      dict created by using find_and_extract_dly
    path_output: Union[Path, str]

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
    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Mean Temp (°C)"] + 273.15, axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "tas"
    da.attrs["standard_name"] = "air_temperature"
    da.attrs["long_name"] = "Near-Surface Air Temperature"
    da.attrs["units"] = "K"
    da.attrs[
        "grid_mapping"
    ] = "regular_lon_lat"  # TODO: There is probably a better CF Convention for point-based data
    da.attrs["comments"] = "station data converted from Mean Temp (°C)"

    # for the first variable, we simply create a dataset from it
    ds = da.to_dataset()

    # import the other variables
    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Max Temp (°C)"] + 273.15, axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "tasmax"
    da.attrs["standard_name"] = "air_temperature maximum"
    da.attrs["long_name"] = "Daily Maximum Near-Surface Temperature maximum"
    da.attrs["units"] = "K"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    da.attrs["comments"] = "station data converted from Max Temp (°C)"
    ds["tasmax"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Min Temp (°C)"] + 273.15, axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "tasmin"
    da.attrs["standard_name"] = "air_temperature minimum"
    da.attrs["long_name"] = "Daily Maximum Near-Surface Temperature minimum"
    da.attrs["units"] = "K"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    da.attrs["comments"] = "station data converted from Min Temp (°C)"
    ds["tasmin"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Heat Deg Days (°C)"], axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "hdd"
    da.attrs["standard_name"] = "heating_degree_days"
    da.attrs[
        "long_name"
    ] = "Number of Degrees Celsius Under a Mean Temperature of 18 °C"
    da.attrs["units"] = "C"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    ds["hdd"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Cool Deg Days (°C)"], axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "cdd"
    da.attrs["standard_name"] = "cooling_degree_days"
    da.attrs["long_name"] = "Number of Degrees Celsius Over a Mean Temperature of 18 °C"
    da.attrs["units"] = "C"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    ds["cdd"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Total Rain (mm)"] / 86400, axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "prlp"
    da.attrs["standard_name"] = "rainfall_flux"
    da.attrs["long_name"] = "Liquid Precipitation"
    da.attrs["units"] = "kg m-2 s-1"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    da.attrs[
        "comments"
    ] = "station data converted from Total Rain (mm) using a density of 1000 kg/m³"
    ds["prlp"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Total Snow (cm)"] / 86400, axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "prsn"
    da.attrs["standard_name"] = "snowfall_flux"
    da.attrs["long_name"] = "Snowfall Flux"
    da.attrs["units"] = "kg m-2 s-1"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    da.attrs[
        "comments"
    ] = "station data converted from Total Snow (cm) using a density of 100 kg/m³"
    ds["prsn"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Total Precip (mm)"] / 86400, axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "pr"
    da.attrs["standard_name"] = "precipitation_flux"
    da.attrs["long_name"] = "Precipitation"
    da.attrs["units"] = "kg m-2 s-1"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    da.attrs[
        "comments"
    ] = "station data converted from Total Precip (mm) using a density of 1000 kg/m³"
    ds["pr"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Snow on Grnd (cm)"] / 100, axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "snd"
    da.attrs["standard_name"] = "surface_snow_thickness"
    da.attrs["long_name"] = "Snow Depth"
    da.attrs["units"] = "m"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    da.attrs["comments"] = "station data converted from Snow on Grnd (cm)"
    ds["snd"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Dir of Max Gust (10s deg)"], axis=1), axis=2
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "sfcWindmax_dir"
    da.attrs["standard_name"] = "wind_gust_from_direction maximum"
    da.attrs[
        "long_name"
    ] = "Direction from which the Daily Maximum Near-Surface Gust Wind Speed maximum Blows"
    da.attrs["units"] = "degree"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    da.attrs["comments"] = "station data converted from Dir of Max Gust (10s deg)"
    ds["sfcWindmax_dir"] = da

    da = xr.DataArray(
        np.expand_dims(
            np.expand_dims(station["data"]["Spd of Max Gust (km/h)"] / 3.6, axis=1),
            axis=2,
        ),
        [
            ("time", time),
            ("lat", [station["latitude"]]),
            ("lon", [station["longitude"]]),
        ],
    )
    da.name = "sfcWindmax"
    da.attrs["standard_name"] = "wind_speed_of_gust maximum"
    da.attrs["long_name"] = "Daily Maximum Near-Surface Gust Wind Speed maximum"
    da.attrs["units"] = "m s-1"
    da.attrs["grid_mapping"] = "regular_lon_lat"
    da.attrs["comments"] = "station data converted from Spd of Max Gust (km/h)"
    ds["sfcWindmax"] = da

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


# This calls _read_single_eccc_dly and appends the data in a single Dict
def _read_multiple_eccc_dly(
    files: Union[List[Union[str, Path]], Generator[Path, None, None]],
    rm_flags: bool = False,
) -> dict:
    """

    Parameters
    ----------
    files : List[Union[str, Path]]
      A list of all the files to append.
    rm_flags : bool
      Removes all the 'Flag' and 'Quality' columns of the ECCC files. Default: False.

    Returns
    -------
    dict
    """

    # Extract the data for each files
    station_meta = None
    datafull = None

    if isinstance(files, GeneratorType):
        files = [f for f in files]

    for i, f in enumerate(files):
        station_meta, data = _read_single_eccc_dly(f)
        if i == 0:
            datafull = data
        else:
            datafull = datafull.append(data, ignore_index=True)

    # change the Date/Time column to a datetime64 type
    datafull["Date/Time"] = pd.to_datetime(datafull["Date/Time"])

    # if wanted, remove the quality and flag columns
    if rm_flags:
        index_quality = [
            i for i, s in enumerate(datafull.columns.values) if "Quality" in s
        ]
        datafull = datafull.drop(datafull.columns.values[index_quality], axis="columns")
        index_flag = [i for i, s in enumerate(datafull.columns.values) if "Flag" in s]
        datafull = datafull.drop(datafull.columns.values[index_flag], axis="columns")

    # combine everything in a single Dict
    station = station_meta
    station["data"] = datafull

    return station


# This is the script that actually reads the CSV files.
# The metadata are saved in a Dict, while the data is returned as a pandas Dataframe.
def _read_single_eccc_dly(file: Union[Path, str]) -> Tuple[dict, pd.DataFrame]:
    """

    Parameters
    ----------
    file : Union[Path, str]

    Returns
    -------
    Tuple[dict, pd.DataFrame]
    """
    # Read the whole file
    with open(file, "r") as fi:
        lines = fi.readlines()

    # Find each elements in the header
    search_header = [0] * 9
    search_header[0] = [i for i, s in enumerate(lines) if "Station Name" in s][0]
    search_header[1] = [i for i, s in enumerate(lines) if "Province" in s][0]
    search_header[2] = [i for i, s in enumerate(lines) if "Latitude" in s][0]
    search_header[3] = [i for i, s in enumerate(lines) if "Longitude" in s][0]
    search_header[4] = [i for i, s in enumerate(lines) if "Elevation" in s][0]
    search_header[5] = [i for i, s in enumerate(lines) if "Climate Identifier" in s][0]
    search_header[6] = [i for i, s in enumerate(lines) if "WMO Identifier" in s][0]
    search_header[7] = [i for i, s in enumerate(lines) if "TC Identifier" in s][0]
    search_header[8] = [i for i, s in enumerate(lines) if "Date/Time" in s][
        0
    ]  # This is where the data actually starts

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
            "Data for file {} is not starting on January 1st. Make sure this is what you want!".format(
                file.name
            )
        )

    return station_meta, data
