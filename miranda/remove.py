import logging
from datetime import date
from getpass import getpass
from logging import config
from pathlib import Path
from types import GeneratorType
from typing import List, Optional, Union

import fabric

from .scripting import LOGGING_CONFIG
from .storage import report_file_size
from .utils import creation_date, ingest

config.dictConfig(LOGGING_CONFIG)


def file_emptier(*, file_list: List[Union[str, Path]]) -> None:
    """
    Provided a list of file paths, will open and overwrite them in order to delete data while preserving the file name.

    Parameters
    ----------
    file_list: List[Union[str, Path]]
      List of files to be overwritten

    Returns
    -------
    None
    """

    file_list = ingest(file_list)

    logging.info(
        "Found {} files totalling {}".format(
            len(file_list), report_file_size(file_list)
        )
    )

    for file in file_list:
        logging.warning(f"Overwriting {file}")
        open(file, "w").close()


def delete_by_date(
    *,
    source: Union[str, Path],
    year: Optional[int] = None,
    month: Optional[int] = None,
    day: Optional[int] = None,
    pattern: Optional[str] = None,
    server: Optional[Union[str, Path]] = None,
    user: Optional[str] = None,
    password: Optional[str] = None,
    date_object: Optional[date] = None,
) -> None:
    """

    Parameters
    ----------
    source: Union[str, Path]
    year: Optional[int]
    month: Optional[int]
    day: Optional[int]
    pattern: Optional[str]
    server: Optional[Union[str, Path]]
    user: Optional[str]
    password: Optional[str]
    date_object: Optional[date]

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

    logging.info(f"Found {len(nc_files)} files totalling {report_file_size(nc_files)}")

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
        "Removed {} files totalling {}".format(
            deleted_files, report_file_size(freed_space)
        )
    )

    if server:
        context.close()

    return


def delete_duplicates(
    *,
    source: Union[str, Path],
    target: Union[str, Path],
    server: Optional[Union[str, Path]],
    user: str = None,
    password: str = None,
    pattern: str = None,
    delete_target_duplicates: bool = False,
) -> None:
    """

    Parameters
    ----------
    source : Union[str, Path]
    target : Union[str, Path]
    server : Optional[Union[str, Path]]
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
        "Found {} files totalling {}".format(
            len(nc_file_duplicates), report_file_size(nc_file_duplicates)
        )
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
        "Removed {} files totalling {}".format(
            deleted_files, report_file_size(freed_space)
        )
    )
    return


def delete_by_variable(
    *,
    target: Union[str, Path, List[Union[str, Path]], GeneratorType] = None,
    variables: List[str] = None,
    server: Optional[str or Path],
    user: Optional[str] = None,
    password: Optional[str] = None,
    file_suffix: Optional[str] = None,
    delete: bool = False,
) -> None:
    """
    Given target location(s), a list of variables and a server address, perform a glob search
     and delete file names starting with the variables identified

    Parameters
    ----------
    target : Union[str, Path, List[Union[str, Path]], GeneratorType]
    variables : List[str]
    server :Optional[Union[str, Path]]
    user : Optional[str]
    password : Optional[str]
    file_suffix : Optional[str]
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
            "Found {} files totalling {}".format(
                len(nc_files), report_file_size(nc_files)
            )
        )

        with connection as context:
            for file in nc_files:
                freed_space += Path(file).stat().st_size
                deleted_files += 1
                if delete:
                    logging.info(f"Deleting file {file.stem}")
                    context.remove(file)

    logging.info(
        "Removed {} files totalling {}".format(
            deleted_files, report_file_size(freed_space)
        )
    )
    return
