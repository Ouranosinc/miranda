import logging
import os
import platform
import sys
from contextlib import contextmanager
from datetime import date
from datetime import datetime as dt
from pathlib import Path
from types import GeneratorType
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional
from typing import Sequence
from typing import Union

KiB = int(pow(2, 10))
MiB = int(pow(2, 20))
GiB = int(pow(2, 30))


def _ingest(files: Union[GeneratorType, List]) -> List:
    if isinstance(files, GeneratorType):
        files = [f for f in files]
    files.sort()
    return files


def creation_date(path_to_file: Union[Path, str]) -> Union[float, date]:
    """
    Try to get the date that a file was created, falling back to when it was last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.

    Parameters
    ----------
    path_to_file: Union[Path, str]

    Returns
    -------
    Union[float, date]
    """
    if platform.system() == "Windows":
        return Path(path_to_file).stat().st_ctime
    else:
        stat = Path(path_to_file).stat()
        try:
            return date.fromtimestamp(stat.st_ctime)
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return date.fromtimestamp(stat.st_mtime)


def read_privileges(location: Union[Path, str]) -> bool:
    """
    Determine whether a user has read privileges to a specific file

    Parameters
    ----------
    location: Union[Path, str]

    Returns
    -------
    bool
      Whether or not the current user shell has read privileges
    """
    if (2, 7) < sys.version_info < (3, 6):
        location = str(location)

    try:
        if Path(location).exists():
            if os.access(location, os.R_OK):
                logging.info(
                    "{} is read OK!".format(dt.now().strftime("%Y-%m-%d %X"), location)
                )
                return True
            else:
                msg = "Ensure read privileges."
                logging.error(msg)
                return False
        else:
            logging.error("{} is an invalid path.")
            return False
    except OSError:
        msg = "Ensure read privileges."
        logging.exception(msg)
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
    return


def find_filepaths(
    source: Union[Path, str, GeneratorType, List[Union[Path, str]]],
    recursive: bool = True,
    file_suffixes: Optional[Union[str, List[str]]] = None,
    **_
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
        file_suffixes = list().append(["*", ".*"])
    elif isinstance(file_suffixes, str):
        file_suffixes = [file_suffixes]

    found = list()
    if isinstance(source, (Path, str)):
        source = [source]

    for location in source:
        for pattern in file_suffixes:
            if "*" not in pattern:
                pattern = "*{}".format(pattern)
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


def make_local_dirs(pathway: Union[str, Path], mode: Union[int, bytes] = 0o777) -> None:
    """Create directories recursively, unless they already exist.

    Parameters
    ----------
    pathway : Union[Path, str]
      Path of folders to create.
    mode : Union[int, bytes]

    Returns
    -------
    None
    """

    pathway = Path(pathway)

    if not pathway.exists():
        try:
            pathway.mkdir(parents=True, mode=mode)
        except OSError:
            raise


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
    elif user_input.lower() == "n":
        return False
    else:
        raise ValueError("{} not in (y, n)".format(user_input))


def verbose_fn(message: str, verbose=True) -> None:
    """Trigger verbose mode.
    Parameters
    ----------
    message : str
    verbose : bool
        flag for whether of not to output the message (default: True).
    """

    if verbose:
        print(message)


def list_paths_with_elements(base_paths: List[str], elements: List[str]) -> List[Dict]:
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
    if not hasattr(base_paths, "__iter__"):
        base_paths = [base_paths]
    base_paths = map(os.path.abspath, base_paths)
    # If elements list is empty, return empty list (end of recursion).
    if not elements:
        return []
    #
    paths_elements = []
    for base_path in base_paths:
        try:
            path_content = os.listdir(base_path)
        except NotADirectoryError:
            continue
        path_content.sort()
        next_base_paths = []
        for path_item in path_content:
            next_base_paths.append(os.path.join(base_path, path_item))
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


def eccc_hourly_variable_metadata(variable_name: str) -> dict:
    """
    Un fonction qui retourne differentes informations en fonction de la variable voulue

    Parameters
    ----------
    variable_name: str

    Returns
    -------
    dict
    """

    if variable_name == "wind_speed":
        var_code = 76
        unites = "km h-1"
        fact_mlt = 1
    elif variable_name == "station_pressure":
        var_code = 77
        unites = "kPa"
        fact_mlt = 0.01
    elif variable_name == "dry_bulb_temperature":
        var_code = 78
        unites = "degC"
        fact_mlt = 0.1
    elif variable_name == "relative_humidity":
        var_code = 80
        unites = "%"
        fact_mlt = 1
    elif variable_name == "freezing_rain":
        var_code = 89
        unites = "1"
        fact_mlt = 1
    elif variable_name == "ice_pellet_presence":
        var_code = 94
        unites = "1"
        fact_mlt = 1
    elif variable_name == "hourly_rainfall":
        var_code = 123
        unites = "mm"
        fact_mlt = 0.1
    elif variable_name == "precipitation_amount":
        var_code = 262
        unites = "mm"
        fact_mlt = 0.1
    else:
        msg = 'Variable name "{}" not recognized'.format(variable_name)
        logging.error(msg)
        raise RuntimeError(msg)

    fact_add = 0.0
    missing_flags = ["M"]
    least_sig_digit = None

    return dict(
        code_var=var_code,
        unites=unites,
        fact_mlt=fact_mlt,
        fact_add=fact_add,
        flag_manquants=missing_flags,
        least_significant_digit=least_sig_digit,
    )


def eccc_cd_hourly_metadata(variable_code: Union[int, str]) -> dict:
    """

    Parameters
    ----------
    variable_code: Union[int, str]

    Returns
    -------
    dict
    """
    ec_hourly_variables = {
        "076": {
            "nc_units": "m s-1",
            "scale_factor": 0.277777778,
            "add_offset": 0,
            "long_name": "Wind Speed (U2A)",
            "standard_name": "wind_speed",
            "nc_name": "windspeed",
        },
        "077": {
            "nc_units": "Pa",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Station Pressure",
            "standard_name": "pressure",
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
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Rainfall",
            "standard_name": "rainfall_amount",
            "nc_name": "rainfall",
        },
        "089": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Freezing Rain",
            "standard_name": "freezing_rain",
            "nc_name": "freezerain",
        },
        "094": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Ice Pellets",
            "standard_name": "ice_pellet_presence",
            "nc_name": "icepellets",
        },
        "123": {
            "nc_units": "%",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Relative Humidity",
            "standard_name": "relative_humidity",
            "nc_name": "hur",
        },
        "262": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 00-60)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_hourly",
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
            "nc_name": "precip_weight_q1",
        },
        "268": {
            "nc_units": "mk m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 30)",
            "standard_name": "precipitation_amount",
            "nc_name": "precip_weight_q2",
        },
        "269": {
            "nc_units": "m",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 45)",
            "standard_name": "precipitation_amount",
            "nc_name": "precip_weight_q3",
        },
        "270": {
            "nc_units": "m",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 60)",
            "standard_name": "precipitation_amount",
            "nc_name": "precip_weight_q4",
        },
        "271": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 00-15)",
            "standard_name": "wind_speed",
            "nc_name": "windspeed_q1",
        },
        "272": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 15-30)",
            "standard_name": "wind_speed",
            "nc_name": "windpseed_q2",
        },
        "273": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 30-45)",
            "standard_name": "wind_speed",
            "nc_name": "windspeed_q3",
        },
        "274": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 45-60)",
            "standard_name": "wind_speed",
            "nc_name": "windspeed_q4",
        },
        "275": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 60)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd",
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
            "nc_name": "snd_q4",
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
            "nc_name": "wnd",
        },
    }
    code = str(variable_code).zfill(3)
    return ec_hourly_variables[code]


def eccc_cf_daily_metadata(variable_code: Union[int, str]) -> dict:
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
            "standard_name": "air_temperature",
            "nc_name": "tasmax",
        },
        "002": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Daily Minimum Temperature",
            "standard_name": "air_temperature",
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
            "long_name": "Total Rainfall",
            "standard_name": "rainfall_accumulation",
            "nc_name": "rainfall",
        },
        "011": {
            "nc_units": "cm",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Snowfall",
            "standard_name": "snowfall_accumulation",
            "nc_name": "snowfall",
        },
        "012": {
            "nc_units": "mm",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation",
            "standard_name": "precipitation_accumulation",
            "nc_name": "precipitation",
        },
        "013": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow on the Ground",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd",
        },
    }
    code = str(variable_code).zfill(3)
    return ec_daily_variables[code]
