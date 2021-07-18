import logging.config
import os
import sys
from contextlib import contextmanager
from datetime import date
from datetime import datetime as dt
from pathlib import Path
from types import GeneratorType
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from . import scripting

KiB = int(pow(2, 10))
MiB = int(pow(2, 20))
GiB = int(pow(2, 30))

logging.config.dictConfig(scripting.LOGGING_CONFIG)

__all__ = [
    "creation_date",
    "eccc_ahccd_metadata",
    "eccc_cf_daily_metadata",
    "eccc_cf_hourly_metadata",
    "find_filepaths",
    "GiB",
    "ingest",
    "KiB",
    "list_paths_with_elements",
    "MiB",
    "read_privileges",
    "set_comparisons",
    "single_item_list",
    "working_directory",
    "yesno_prompt",
]


def ingest(files: Union[GeneratorType, List]) -> List:
    if isinstance(files, GeneratorType):
        files = [f for f in files]
    files.sort()
    return files


def creation_date(path_to_file: Union[Path, str]) -> Union[float, date]:
    """
    Try to get the date that a file was created, falling back to when it was last modified if that isn't possible.
    See https://stackoverflow.com/a/39501288/1709587 for explanation.

    Parameters
    ----------
    path_to_file: Union[Path, str]

    Returns
    -------
    Union[float, date]
    """
    if os.name == "nt":
        return Path(path_to_file).stat().st_ctime

    stat = Path(path_to_file).stat()
    try:
        return date.fromtimestamp(stat.st_ctime)
    except AttributeError:
        # We're probably on Linux. No easy way to get creation dates here,
        # so we'll settle for when its content was last modified.
        return date.fromtimestamp(stat.st_mtime)


def read_privileges(location: Union[Path, str], strict: bool = False) -> bool:
    """
    Determine whether a user has read privileges to a specific file

    Parameters
    ----------
    location: Union[Path, str]
    strict: bool

    Returns
    -------
    bool
      Whether or not the current user shell has read privileges
    """
    if (2, 7) < sys.version_info < (3, 6):
        location = str(location)

    msg = str()
    try:
        if Path(location).exists():
            if os.access(location, os.R_OK):
                msg = "{} is read OK!".format(location)
                logging.info(msg)
                return True
            msg = "Ensure read privileges for `{}`.".format(location)
        else:
            msg = "`{}` is an invalid path.".format(location)
        raise OSError

    except OSError:
        logging.exception(msg)
        if strict:
            raise
        return False


@contextmanager
def working_directory(directory: Union[str, Path]) -> None:
    """
    This function momentarily changes the working directory within the context and reverts to the file working directory
    when the code block it is acting upon exits

    Parameters
    ----------
    directory: Union[str, Path]

    Returns
    -------
    None

    """
    owd = os.getcwd()

    if (2, 7) < sys.version_info < (3, 6):
        directory = str(directory)

    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)


def find_filepaths(
    source: Union[Path, str, GeneratorType, List[Union[Path, str]]],
    recursive: bool = True,
    file_suffixes: Optional[Union[str, List[str]]] = None,
    **_,
) -> List[Path]:
    """

    Parameters
    ----------
    source : Union[Path, str, GeneratorType, List[Union[Path, str]]]
    recursive : bool
    file_suffixes: List[str]

    Returns
    -------
    List[Path]
    """

    if file_suffixes is None:
        file_suffixes = ["*", ".*"]
    elif isinstance(file_suffixes, str):
        file_suffixes = [file_suffixes]

    found = list()
    if isinstance(source, (Path, str)):
        source = [source]

    for location in source:
        for pattern in file_suffixes:
            if "*" not in pattern:
                pattern = "*{}*".format(pattern)
            if recursive:
                found.extend([f for f in Path(location).expanduser().rglob(pattern)])
            elif not recursive:
                found.extend([f for f in Path(location).expanduser().glob(pattern)])
            else:
                raise ValueError("Recursive: {}".format(recursive))

    if (2, 7) < sys.version_info < (3, 6):
        found = [str(f) for f in found]

    return found


def single_item_list(iterable: Iterable) -> bool:
    """
    See: https://stackoverflow.com/a/16801605/7322852

    Parameters
    ----------
    iterable: Iterable

    Returns
    -------
    bool

    """
    iterator = iter(iterable)
    has_true = any(iterator)  # consume from "i" until first true or it's exhausted
    has_another_true = any(
        iterator
    )  # carry on consuming until another true value / exhausted
    return has_true and not has_another_true  # True if exactly one true found


def set_comparisons(set1: Sequence, set2: Sequence) -> bool:
    """Compare two sequences of non hashable objects as if they were sets.

    Parameters
    ----------
    set1 : Sequence
      First sequence of objects.
    set2 : Sequence
      Second sequence of objects.

    Returns
    -------
    out : bool
      True if two sets are identical, that is, they contain the same elements.
    """

    for item1 in set1:
        if item1 not in set2:
            return False
    for item2 in set2:
        if item2 not in set1:
            return False
    return True


########################################################################################


def yesno_prompt(query: str) -> bool:
    """Prompt user for a yes/no answer.
    Parameters
    ----------
    query : str
        the yes/no question to ask the user.

    Returns
    -------
    out : bool
        True (yes) or False (otherwise).
    """

    user_input = input("{} (y/n) ".format(query))
    if user_input.lower() == "y":
        return True
    if user_input.lower() == "n":
        return False
    raise ValueError("{} not in (y, n)".format(user_input))


def list_paths_with_elements(
    base_paths: Union[str, List[str]], elements: List[str]
) -> List[Dict]:
    """List a given path structure.
    Parameters
    ----------
    base_paths : List[str]
        list of paths from which to start the search.
    elements : List[str]
        ordered list of the expected elements.
    Returns
    -------
    out : List[Dict]
        the keys are 'path' and each of the members of the given elements,
        the path is the absolute path.
    Notes
    -----
    Suppose you have the following structure:
    /base_path/{color}/{shape}
    The resulting list would look like:
    [{'path':/base_path/red/square, 'color':'red', 'shape':'square'},
     {'path':/base_path/red/circle, 'color':'red', 'shape':'circle'},
     {'path':/base_path/blue/triangle, 'color':'blue', 'shape':'triangle'},
     ...
    ]
    Obviously, 'path' should not be in the input list of elements.
    """

    # Make sure the base_paths input is a list of absolute path
    paths = list()
    if not hasattr(base_paths, "__iter__"):
        paths.append(base_paths)
    paths = map(os.path.abspath, base_paths)
    # If elements list is empty, return empty list (end of recursion).
    if not elements:
        return list()

    paths_elements = list()
    for base_path in paths:
        try:
            path_content = [f for f in Path(base_path).iterdir()]
        except NotADirectoryError:
            continue
        path_content.sort()
        next_base_paths = []
        for path_item in path_content:
            next_base_paths.append(base_path.joinpath(path_item))
        next_pe = list_paths_with_elements(next_base_paths, elements[1:])
        if next_pe:
            for i, one_pe in enumerate(next_pe):
                relative_path = next_pe[i]["path"].replace(base_path, "", 1)
                new_element = relative_path.split("/")[1]
                next_pe[i][elements[0]] = new_element
            paths_elements.extend(next_pe)
        elif len(elements) == 1:
            for my_path, my_item in zip(next_base_paths, path_content):
                paths_elements.append({"path": my_path, elements[0]: my_item})
    return paths_elements


def eccc_cf_hourly_metadata(
    variable_code: Union[int, str]
) -> Dict[str, Union[int, float]]:
    """

    Parameters
    ----------
    variable_code: Union[int, str]

    Returns
    -------
    dict
    """
    ec_hourly_variables = {
        "061": {
            "nc_units": "W s m-2",  # FIXME: Original units are MJ m-2 h-1, how best to convert?
            "scale_factor": None,  # Not sure how best to proceed here.
            "add_offset": 0,
            "long_name": "RF1 Global Solar Radiation",
            "standard_name": "solar_radiation_flux",
            "nc_name": "rf1_radiation",
        },
        "071": {
            "nc_units": "m",
            "scale_factor": 30,
            "add_offset": 0,
            "long_name": "Ceiling height of lowest layer of clouds",
            "standard_name": "ceiling_cloud_height",
            "nc_name": "ceiling_hgt",
        },
        "072": {
            "nc_units": "m",
            "scale_factor": 100,
            "add_offset": 0,
            "long_name": "Visibility",
            "standard_name": "visibility_in_air",
            "nc_name": "visibility",
        },
        "073": {
            "nc_units": "Pa",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Sea Level Pressure",
            "standard_name": "air_pressure_at_mean_sea_level",
            "nc_name": "psl",
        },
        "074": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Dew Point Temperature",
            "standard_name": "dew_point_temperature",
            "nc_name": "tds",
        },
        "075": {
            "nc_units": "degree",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Wind Direction at 2 m (U2A Anemometer) (16 pts)",
            "standard_name": "wind_direction_u2a",
            "nc_name": "wind_dir_u2a_16",
        },
        "076": {
            "nc_units": "m s-1",
            "scale_factor": 0.277777778,
            "add_offset": 0,
            "long_name": "Wind Speed (U2A Anemometer)",
            "standard_name": "wind_speed_u2a",
            "nc_name": "wind_speed_u2a",
        },
        "077": {
            "nc_units": "Pa",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Station Pressure",
            "standard_name": "atmospheric_pressure",
            "nc_name": "pressure",
        },
        "078": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Dry Bulb Temperature",
            "standard_name": "dry_bulb_temperature",
            "nc_name": "tas_dry",
        },
        "079": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Wet Bulb temperature",
            "standard_name": "wet_bulb_temperature",
            "nc_name": "tas_wet",
        },
        "080": {
            "nc_units": "%",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Relative Humidity",
            "standard_name": "relative_humidity",
            "nc_name": "hur",
        },
        "089": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Freezing Rain",
            "standard_name": "freezing_rain",
            "nc_name": "freeze_rain",
        },
        "094": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Ice Pellets",
            "standard_name": "ice_pellet_presence",
            "nc_name": "ice_pellets",
        },
        "107": {
            "nc_units": "%",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Lowest cloud layer opacity",
            "standard_name": "low_type_cloud_opacity_fraction",
            "nc_name": "cloud_opac",
        },
        "108": {
            "nc_units": "%",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Lowest cloud layer amount or condition",
            "standard_name": "low_type_cloud_area_fraction",
            "nc_name": "cloud_frac",
        },
        "109": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Lowest cloud layer type",
            "standard_name": "low_type_cloud_type",
            "nc_name": "low_cloud_type",
        },
        "110": {
            "nc_units": "m",
            "scale_factor": 30,
            "add_offset": 0,
            "long_name": "Lowest cloud layer height",
            "standard_name": "low_type_cloud_height",
            "nc_name": "low_cloud_hgt",
        },
        "123": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Rainfall",
            "standard_name": "rainfall_amount",
            "nc_name": "rainfall",
        },
        "133": {
            "nc_units": "s",
            "scale_factor": 3600,
            "add_offset": 0,
            "long_name": "Sunshine",
            "standard_name": "duration_of_sunshine",
            "nc_name": "sun",
        },
        "156": {
            "nc_units": "degree",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Wind Direction at 2 m (U2A Anemometer) (36 pts)",
            "standard_name": "wind_direction_u2a",
            "nc_name": "wind_dir_u2a_36",
        },
        "262": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 00-60)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation",
        },
        "263": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 00-15)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_q1",
        },
        "264": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 15-30)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_q2",
        },
        "265": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 30-45)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_q3",
        },
        "266": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 45-60)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_q4",
        },
        "267": {
            "nc_units": "kg m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 15)",
            "standard_name": "precipitation_amount",
            "nc_name": "precipitation_weight_q1",
        },
        "268": {
            "nc_units": "kg m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 30)",
            "standard_name": "precipitation_amount",
            "nc_name": "precipitation_weight_q2",
        },
        "269": {
            "nc_units": "kg m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 45)",
            "standard_name": "precipitation_amount",
            "nc_name": "precipitation_weight_q3",
        },
        "270": {
            "nc_units": "kg m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 60)",
            "standard_name": "precipitation_amount",
            "nc_name": "precipitation_weight_q4",
        },
        "271": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 00-15)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed_q1",
        },
        "272": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 15-30)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed_q2",
        },
        "273": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 30-45)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed_q3",
        },
        "274": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 45-60)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed_q4",
        },
        "275": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 60)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd_q4",
        },
        "276": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 15)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd_q1",
        },
        "277": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 30)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd_q2",
        },
        "278": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 45)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd_q3",
        },
        "279": {
            "nc_units": "degree",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Wind Direction at 2 m (minutes 50-60)",
            "standard_name": "wind_direction",
            "nc_name": "wind_dir",
        },
        "280": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 50-60)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed",
        },
    }
    code = str(variable_code).zfill(3)
    if code in ["061"]:
        raise NotImplementedError()
    try:
        variable = ec_hourly_variables[code]
        variable["missing_flags"] = "M"
        variable["least_significant_digit"] = ""
    except KeyError:
        logging.error("Hourly variable `{}` not supported.".format(code))
        raise
    return variable


def eccc_cf_daily_metadata(
    variable_code: Union[int, str]
) -> Dict[str, Union[int, float]]:
    """

    Parameters
    ----------
    variable_code: Union[int, str]

    Returns
    -------
    dict
    """
    ec_daily_variables = {
        "001": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Daily Maximum Temperature",
            "standard_name": "air_temperature_maximum",
            "nc_name": "tasmax",
        },
        "002": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Daily Minimum Temperature",
            "standard_name": "air_temperature_minimum",
            "nc_name": "tasmin",
        },
        "003": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Daily Mean Temperature",
            "standard_name": "air_temperature",
            "nc_name": "tas",
        },
        "010": {
            "nc_units": "mm",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Daily Total Rainfall",
            "standard_name": "liquid_precipitation_flux",
            "nc_name": "prlptot",
        },
        "011": {
            "nc_units": "cm",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Daily Total Snowfall",
            "standard_name": "solid_precipitation_flux",
            "nc_name": "prsntot",
        },
        "012": {
            "nc_units": "mm",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Daily Total Precipitation",
            "standard_name": "precipitation_flux",
            "nc_name": "prcptot",
        },
        "013": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow on the Ground",
            "standard_name": "surface_snow_thickness",
            "nc_name": "sndtot",
        },
        "014": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Thunderstorms",
            "standard_name": "thunderstorm_presence",
            "nc_name": "thunder",
        },
        "015": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Freezing rain or drizzle",
            "standard_name": "freeze_rain_drizzle_presence",
            "nc_name": "freezing_rain_drizzle",
        },
        "016": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Hail",
            "standard_name": "hail_presence",
            "nc_name": "hail",
        },
        "017": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Fog or Ice Fog",
            "standard_name": "fog_ice_fog_presence",
            "nc_name": "fog_ice_fog",
        },
        "018": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Smoke or Haze",
            "standard_name": "smoke_haze_presence",
            "nc_name": "smoke_haze",
        },
        "019": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Blowing Dust or Sand",
            "standard_name": "blowing_dust_sand_presence",
            "nc_name": "blowing_dust_sand",
        },
        "020": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Blowing snow",
            "standard_name": "blowing_snow_presence",
            "nc_name": "blow_snow",
        },
        "021": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Wind speed >= 28 Knots",
            "standard_name": "wind_exceeding_28_knots",
            "nc_name": "wind_gt_28kt",
        },
        "022": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Wind speed >= 34 Knots",
            "standard_name": "wind_exceeding_34_knots",
            "nc_name": "wind_gt_34kt",
        },
        "023": {
            "nc_units": "degree",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Direction of extreme gust (16 pts) to December 1976",
            "standard_name": "wind_to_direction",
            "nc_name": "gust_dir",
        },
        "024": {
            "nc_units": "m s-1",
            "scale_factor": 0.2777778,
            "add_offset": 0,
            "long_name": "Speed of extreme gust",
            "standard_name": "wind_speed_of_gust",
            "nc_name": "gust_speed",
        },
        "025": {
            "nc_units": "h",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "UTC hour of extreme gust",
            "standard_name": "hour_of_extreme_gust",
            "nc_name": "gust_hour",
        },
    }
    code = str(variable_code).zfill(3)
    try:
        variable = ec_daily_variables[code]
        variable["missing_flags"] = "M"
        variable["least_significant_digit"] = ""
    except KeyError:
        logging.error("Daily variable `{}` not supported.".format(code))
        raise
    return variable


def eccc_ahccd_metadata(
    code: str, gen: int
) -> (dict, dict, List[Tuple[int, int]], int):
    """

    Parameters
    ----------
    code: str
    gen: int

    Returns
    -------
    dict, dict, List[Tuple[(int, int)]], int
    """

    ec_ahccd_attrs = dict(
        dx=dict(
            variable="tasmax",
            units="deg C",
            standard_name="air_temperature",
            long_name="Near-Surface Maximum Daily Air Temperature",
            comment=f"ECCC {gen} Generation of Adjusted and Homogenized Temperature Data",
        ),
        dn=dict(
            variable="tasmin",
            units="deg C",
            standard_name="air_temperature",
            long_name="Near-Surface Minimum Daily Air Temperature",
            comment=f"ECCC {gen} Generation of Adjusted and Homogenized Temperature Data",
        ),
        dm=dict(
            variable="tas",
            units="deg C",
            standard_name="air_temperature",
            long_name="Near-Surface Daily Mean Air Temperature",
            comment=f"ECCC {gen} Generation of Adjusted and Homogenized Temperature Data",
        ),
        dt=dict(
            variable="pr",
            units="mm day-1",
            standard_name="precipitation_flux",
            long_name="Daily Total Precipitation",
            comment=f"ECCC {gen} Generation of Adjusted and Homogenized Precipitation Data",
        ),
        ds=dict(
            variable="prsn",
            units="mm day-1",
            standard_name="snowfall_flux",
            long_name="Daily Snowfall",
            comment=f"ECCC {gen} Generation of Adjusted and Homogenized Precipitation Data",
        ),
        dr=dict(
            variable="prlp",
            units="mm day-1",
            standard_name="rainfall_flux",
            long_name="Daily Rainfall",
            comment=f"ECCC {gen} Generation of Adjusted and Homogenized Precipitation Data",
        ),
    )
    try:
        variable = ec_ahccd_attrs[code]
        variable["missing_flags"] = "M"
        if variable["variable"].startswith("tas"):
            variable["NaN_value"] = "-9999.9"
            column_names = [
                "Prov",
                "Station name",
                "stnid",
                "FromYear",
                "FromMonth",
                "ToYear",
                "ToMonth",
                "lat",
                "long",
                "elev",
                "stns joined",
            ]
            column_spaces = [(0, 5), (5, 6), (6, 8), (8, 9)]
            ii = 9
            for i in range(1, 32):
                column_spaces.append((ii, ii + 7))
                ii += 7
                column_spaces.append((ii, ii + 1))
                ii += 1
            header_row = 3

        elif variable["variable"].startswith("pr"):
            variable["NaN_value"] = "-9999.99"
            column_names = [
                "Prov",
                "Station name",
                "stnid",
                "beg yr",
                "beg mon",
                "end yr",
                "end mon",
                "lat (deg)",
                "long (deg)",
                "elev (m)",
                "stns joined",
            ]
            column_spaces = [(0, 4), (4, 5), (5, 7), (7, 8)]
            ii = 8
            for i in range(1, 32):
                column_spaces.append((ii, ii + 8))
                ii += 8
                column_spaces.append((ii, ii + 1))
                ii += 1
            header_row = 0

        else:
            raise KeyError

        column_names = [
            col.lower().split("(")[0].replace("%", "pct_").replace(" ", "_")
            for col in list(column_names)
        ]

        generation = {1: "First", 2: "Second", 3: "Third"}.get(gen)
        if gen == 1:
            _citation = "Vincent, L.A., M.M. Hartwell and X.L. Wang, 2020: A Third Generation of Homogenized "
            "Temperature for Trend Analysis and Monitoring Changes in Canada’s Climate. "
            "Atmosphere-Ocean. https://doi.org/10.1080/07055900.2020.1765728"
        elif gen == 2:
            _citation = "Mekis, É and L.A. Vincent, 2011: An overview of the second generation adjusted daily "
            "precipitation dataset for trend analysis in Canada. Atmosphere-Ocean 49(2), "
            "163-177 doi:10.1080/07055900.2011.583910"
        else:
            msg = f"Generation '{gen}' not supported."
            raise NotImplementedError(msg)

        global_attrs = dict(
            title=f"{generation} Generation of Homogenized Daily {variable} for Canada Update to December 2019",
            history=f"{dt.today().strftime('%Y-%m-%d')}: Convert from original format to NetCDF",
            type="station_obs",
            institute="Environment and Climate Change Canada",
            institute_id="ECCC",
            dataset_id=f"AHCCD_gen{gen}_day_{variable}",
            frequency="day",
            licence_type="permissive",
            licence="https:/open.canada.ca/en/open-government-licence-canada",
            citation=_citation,
        )

    except KeyError as e:
        msg = f"AHCCD variable '{code}' or generation '{gen}' not supported."
        logging.error(msg)
        raise NotImplementedError(msg) from e

    return variable, column_names, column_spaces, header_row, global_attrs
