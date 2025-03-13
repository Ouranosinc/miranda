from __future__ import annotations

import logging.config
import os
from pathlib import Path
from types import GeneratorType

import netCDF4 as nc  # noqa

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = [
    "discover_data",
]


# FIXME: How are these two functions different?
def discover_data(
    input_files: str | os.PathLike | list[str | os.PathLike] | GeneratorType,
    suffix: str = "nc",
    recurse: bool = True,
) -> list[Path] | GeneratorType:
    """Discover data.

    Parameters
    ----------
    input_files : str, pathlib.Path, list of str or Path, or GeneratorType
        Path or string to a file, a folder, or a generator of paths.
    suffix : str
        File-ending suffix to search for. Default: "nc".
    recurse : bool
        Whether to recurse through folders or not. Default: True.

    Returns
    -------
    list of pathlib.Path or GeneratorType of pathlib.Path

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
