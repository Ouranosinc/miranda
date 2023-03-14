import logging.config
import os
import sys
from datetime import date
from pathlib import Path
from types import GeneratorType
from typing import List, Optional, Union

import netCDF4
import zarr

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = [
    "get_chunks_on_disk",
    "discover_data",
    "creation_date",
    "read_privileges",
    "find_filepaths",
]


def get_chunks_on_disk(file: Union[os.PathLike, str]) -> dict:
    """

    Parameters
    ----------
    file : str or os.PathLike
        File to be examined. Supports NetCDF and Zarr.

    Returns
    -------
    dict
    """
    chunks = dict()
    file = Path(file)

    if file.suffix.lower() in [".nc", ".nc4"]:
        with netCDF4.Dataset(file) as ds:
            for v in ds.variables:
                chunks[v] = dict()
                for ii, dim in enumerate(ds[v].dimensions):
                    chunks[v][dim] = ds[v].chunking()[ii]
    elif file.suffix.lower() == "zarr" and file.is_dir():
        with zarr.open(file, "r") as ds:  # noqa
            for v in ds.arrays():
                # Check if variable is chunked
                if v[1]:
                    chunks[v[0]] = v[1]
    else:
        raise NotImplementedError(f"File type: {file.suffix}.")
    return chunks


def discover_data(
    input_files: Union[str, os.PathLike, List[Union[str, os.PathLike]], GeneratorType],
    suffix: str = ".nc",
    recurse: bool = True,
) -> Union[List[os.PathLike], GeneratorType]:
    """Discover data.

    Parameters
    ----------
    input_files: str or Path or List[Union[str, Path]] or GeneratorType
      Path or string to a file, a folder, or a generator of paths.
    suffix: str
      File-ending suffix to search for. Default: ".nc".
    recurse: bool
      Whether to recurse through folders or not. Default: True.

    Returns
    -------
    list or generator of Path

    Warnings
    --------
    Recursion through ".zarr" files is explicitly disabled. Recursive globs and generators will not be expanded/sorted.

    """
    if isinstance(input_files, (Path, str)):
        input_files = Path(input_files)
        if input_files.is_dir():
            if suffix.endswith("zarr") or not recurse:
                input_files = sorted(list(input_files.glob(suffix)))
            else:
                input_files = input_files.rglob(suffix)
    elif isinstance(input_files, list):
        input_files = sorted(Path(p) for p in input_files)
    elif isinstance(input_files, GeneratorType):
        pass
    else:
        raise NotImplementedError(f"input_files: {type(input_files)}")
    return input_files


def creation_date(path_to_file: Union[Path, str]) -> Union[float, date]:
    """Try to get the date that a file was created, falling back to when it was last modified if that isn't possible.

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
    """Determine whether a user has read privileges to a specific file.

    Parameters
    ----------
    location: Union[Path, str]
    strict: bool

    Returns
    -------
    bool
      Whether the current user shell has read privileges
    """
    msg = ""
    try:
        if Path(location).exists():
            if os.access(location, os.R_OK):
                msg = f"{location} is read OK!"
                logging.info(msg)
                return True
            msg = f"Ensure read privileges for `{location}`."
        else:
            msg = f"`{location}` is an invalid path."
        raise OSError()

    except OSError:
        logging.exception(msg)
        if strict:
            raise
        return False


def find_filepaths(
    source: Union[Path, str, GeneratorType, List[Union[Path, str]]],
    recursive: bool = True,
    file_suffixes: Optional[Union[str, List[str]]] = None,
    **_,
) -> List[Path]:
    """Find all available filepaths at a given source.

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
                pattern = f"*{pattern}*"
            if recursive:
                found.extend([f for f in Path(location).expanduser().rglob(pattern)])
            elif not recursive:
                found.extend([f for f in Path(location).expanduser().glob(pattern)])
            else:
                raise ValueError(f"Recursive: {recursive}")

    if (2, 7) < sys.version_info < (3, 6):
        found = [str(f) for f in found]

    return found
