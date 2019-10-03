import logging
from datetime import date
from getpass import getpass
from pathlib import Path
from types import GeneratorType
from typing import List
from typing import Optional

import fabric

from miranda.storage import report_file_size
from miranda.utils import creation_date


def file_emptier(files: List[str or Path]) -> None:
    for f in files:
        logging.warning("Overwriting {}".format(f))
        open(f, "w").close()


def delete_by_date(
    *,
    source: str or Path,
    server: Optional[str or Path],
    username: str,
    password: str,
    year: int,
    month: int,
    day: int,
    pattern: str = None,
    dt_object: Optional[date] = None
) -> None:

    date_selected = date(year, month, day) or dt_object
    glob_pattern = pattern or "*.nc"

    nc_files = Path(source).rglob(glob_pattern)
    nc_files = list(nc_files)
    nc_files.sort()

    context = None
    if server:
        connection = fabric.Connection(
            host=server, user=username, connect_kwargs=dict(password=password)
        )
        context = connection.open()

    freed_space = 0
    deleted_files = 0

    for file in nc_files:
        if creation_date(file) == date_selected:
            freed_space += Path(file).stat().st_size
            logging.info("Deleting {}".format(file.name))
            if server:
                context.remove(file)
            else:
                file.unlink()
            deleted_files += 1
    if server:
        context.close()
    return


def delete_duplicates(
    *,
    source: str or Path = None,
    target: str or Path = None,
    server: Optional[str or Path],
    user: str = None,
    password: str = None,
    pattern: str = None,
    delete_target_duplicates: bool = False
) -> None:

    user = user or input("Username:")
    password = password or getpass("Password:")
    glob_pattern = pattern or "*.nc"

    connection = fabric.Connection(
        host=server, user=user, connect_kwargs=dict(password=password)
    )

    source = source or input("Source files:")
    target = target or input("Target files:")

    nc_files_source = Path(source).rglob(glob_pattern)
    nc_files_source = {f.stem for f in nc_files_source}
    nc_files_target = Path(target).rglob(glob_pattern)

    nc_file_duplicates = []
    for f in nc_files_target:
        if f.name in nc_files_source:
            logging.info("Duplicate found: {}".format(f.name))
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
                logging.info("Deleting {}".format(dup.name))
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
    target: list or GeneratorType = None,
    variables: list = None,
    server: Optional[str or Path],
    user: str = None,
    password: str = None,
    file_suffix: str = None,
    delete=False
) -> None:
    """
    Given a target location, a list of variables and a server address, perform a glob search
     and delete file names starting with the variables identified
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
        nc_files = Path(target).rglob("{}*".format(var, glob_suffix))
        nc_files = list(Path(f) for f in nc_files)
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
                    logging.info("Deleting file {}".format(file.stem))
                    context.remove(file)

    logging.info(
        "Removed {} files totalling {}".format(
            deleted_files, report_file_size(freed_space)
        )
    )
    return
