import logging
from datetime import date
from logging import config
from pathlib import Path
from typing import List, Union

from miranda.io import creation_date, find_filepaths
from miranda.scripting import LOGGING_CONFIG

__all__ = ["select_by_date_modified"]

logging.config.dictConfig(LOGGING_CONFIG)


def select_by_date_modified(
    source: Union[Path, str],
    year: int,
    month: int,
    day: int,
    pattern: str = None,
    datetime: date = None,
) -> List:
    """

    Parameters
    ----------
    source
    year
    month
    day
    pattern
    datetime

    Returns
    -------

    """

    date_selected = date(year, month, day) or datetime

    files = find_filepaths(source, file_suffixes=pattern)

    selected_files = list()
    for file in files:
        if creation_date(file) == date_selected:
            logging.info(f"Selecting {file}.")
            selected_files.append(file)

    return selected_files
