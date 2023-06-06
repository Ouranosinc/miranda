"""Remote File Removal Operations module."""
from __future__ import annotations

import logging.config
import warnings
from datetime import date
from getpass import getpass
from pathlib import Path
from types import GeneratorType

from miranda.io.utils import creation_date
from miranda.scripting import LOGGING_CONFIG
from miranda.storage import report_file_size

logging.config.dictConfig(LOGGING_CONFIG)

try:
    import fabric  # noqa
except ImportError:
    warnings.warn(
        f"{__name__} functions require additional dependencies. "
        "Please install them with `pip install miranda[remote]`."
    )


__all__ = [
    "delete_by_date",
    "delete_by_variable",
    "delete_duplicates",
    "file_emptier",
]


def file_emptier(*, file_list: list[str | Path] | GeneratorType) -> None:
    """Open and overwrite a list of file paths in order to delete data while preserving the file name.

    Parameters
    ----------
    file_list : list of str or Path, or GeneratorType
        List of files to be overwritten

    Returns
    -------
    None
    """
    file_list = sorted([Path(f) for f in file_list])

    logging.info(
        f"Found {len(file_list)} files totalling {report_file_size(file_list)}."
    )

    for file in file_list:
        logging.warning(f"Overwriting {file}")
        open(file, "w").close()


def delete_by_date(
    *,
    source: str | Path,
    year: int | None = None,
    month: int | None = None,
    day: int | None = None,
    pattern: str | None = None,
    server: str | Path | None = None,
    user: str | None = None,
    password: str | None = None,
    date_object: date | None = None,
) -> None:
    """Remove a selection of files based on a given date of last modification.

    Parameters
    ----------
    source: str or Path
    year: int, optional
    month: int, optional
    day: int, optional
    pattern: str, optional
    server: str or Path, optional
    user: str, optional
    password: str, optional
    date_object: date, optional

    Returns
    -------
    None
    """
    user = user or input("Username:")
    password = password or getpass("Password:")

    if year and month and day:
        date_selected = date(year, month, day)
    elif date_object:
        date_selected = date_object
    else:
        raise ValueError

    glob_pattern = pattern or "*.nc"

    nc_files = Path(source).rglob(glob_pattern)
    nc_files = list(nc_files)
    nc_files.sort()

    logging.info(f"Found {len(nc_files)} files totalling {report_file_size(nc_files)}.")

    context = None
    if server:
        connection = fabric.Connection(
            host=server, user=user, connect_kwargs=dict(password=password)
        )
        context = connection.open()

    freed_space = 0
    deleted_files = 0

    for file in nc_files:
        if creation_date(file) == date_selected:
            freed_space += Path(file).stat().st_size
            logging.info(f"Deleting {file.name}")
            if context:
                context.remove(file)
            else:
                file.unlink()
            deleted_files += 1

    logging.info(
        f"Removed {deleted_files} files totalling {report_file_size(freed_space)}"
    )

    if server:
        context.close()

    return


def delete_duplicates(
    *,
    source: str | Path,
    target: str | Path,
    server: str | Path | None = None,
    user: str = None,
    password: str = None,
    pattern: str = None,
    delete_target_duplicates: bool = False,
) -> None:
    """

    Parameters
    ----------
    source : str or Path
    target : str or Path
    server : str or Path, optional
    user: str
    password : str
    pattern: str
    delete_target_duplicates : bool

    Returns
    -------
    None
    """
    user = user or input("Username:")
    password = password or getpass("Password:")
    glob_pattern = pattern or "*.nc"

    connection = fabric.Connection(
        host=server, user=user, connect_kwargs=dict(password=password)
    )

    nc_files_source = Path(source).rglob(glob_pattern)
    nc_files_source = {f.stem for f in nc_files_source}
    nc_files_target = Path(target).rglob(glob_pattern)

    nc_file_duplicates = []
    for f in nc_files_target:
        if f.name in nc_files_source:
            logging.info(f"Duplicate found: {f.name}")
            nc_file_duplicates.append(f)

    nc_file_duplicates.sort()
    logging.info(
        f"Found {len(nc_file_duplicates)} files totalling {report_file_size(nc_file_duplicates)}"
    )

    freed_space = 0
    deleted_files = 0
    if delete_target_duplicates:
        with connection as context:
            for dup in nc_file_duplicates:
                freed_space += Path(dup).stat().st_size
                logging.info(f"Deleting {dup.name}")
                context.remove(dup)
                deleted_files += 1

    logging.info(
        f"Removed { deleted_files} files totalling {report_file_size(freed_space)}."
    )
    return


def delete_by_variable(
    *,
    target: str | Path | list[str | Path] | GeneratorType = None,
    variables: list[str],
    server: str | Path | None = None,
    user: str | None = None,
    password: str | None = None,
    file_suffix: str | None = None,
    delete: bool = False,
) -> None:
    """Delete according to variable name.

    Given target location(s), a list of variables and a server address, perform a glob search
    and delete file names starting with the variables identified

    Parameters
    ----------
    target : str, Path, list of str or Path, or GeneratorType]
    variables : list of str
    server : str or Path, optional
    user : str, optional
    password : str, optional
    file_suffix : str, optional
    delete : bool

    Returns
    -------
    None
    """
    user = user or input("Username:")
    password = password or getpass("Password:")

    connection = fabric.Connection(
        host=server, user=user, connect_kwargs=dict(password=password)
    )

    freed_space = 0
    deleted_files = 0
    for var in variables:
        glob_suffix = file_suffix or ".nc"

        if isinstance(target, (GeneratorType, list)):
            found = list()
            for location in target:
                found.extend([f for f in Path(location).rglob(f"{var}*{glob_suffix}")])
        else:
            found = Path(target).rglob(f"{var}*{glob_suffix}")

        nc_files = [Path(f) for f in found]
        nc_files.sort()

        logging.info(
            f"Found {len(nc_files)} files totalling {report_file_size(nc_files)}"
        )

        with connection as context:
            for file in nc_files:
                freed_space += Path(file).stat().st_size
                deleted_files += 1
                if delete:
                    logging.info(f"Deleting file {file.stem}")
                    context.remove(file)

    logging.info(
        f"Removed {deleted_files} files totalling {report_file_size(freed_space)}"
    )
    return
