import logging.config
import os
from pathlib import Path
from types import GeneratorType
from typing import List, Optional, Union

import netCDF4 as nc  # noqa

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = [
    "discover_data",
    "find_filepaths",
]


# FIXME: How are these two functions different?
def discover_data(
    input_files: Union[str, os.PathLike, List[Union[str, os.PathLike]], GeneratorType],
    suffix: str = "nc",
    recurse: bool = True,
) -> Union[List[Path], GeneratorType]:
    """Discover data.

    Parameters
    ----------
    input_files: str or Path or List[Union[str, Path]] or GeneratorType
      Path or string to a file, a folder, or a generator of paths.
    suffix: str
      File-ending suffix to search for. Default: "nc".
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
                input_files = sorted(list(input_files.glob(f"*.{suffix}")))
            else:
                input_files = input_files.rglob(f"*.{suffix}")
        if input_files.is_file():
            logging.warning(
                "Data discovery yielded a single file. Casting to `list[Path]`."
            )
            input_files = [input_files]
    elif isinstance(input_files, list):
        input_files = sorted(Path(p) for p in input_files)
    elif isinstance(input_files, GeneratorType):
        logging.warning(
            "A Generator was passed to `discover_data`. Passing object along..."
        )
        pass
    else:
        raise NotImplementedError(f"input_files: {type(input_files)}")
    return input_files


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

    return found
