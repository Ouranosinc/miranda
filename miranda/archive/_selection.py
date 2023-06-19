from __future__ import annotations

import logging
from datetime import datetime
from logging import config
from pathlib import Path

from miranda.io import find_filepaths
from miranda.io.utils import creation_date
from miranda.scripting import LOGGING_CONFIG

__all__ = ["select_by_date_modified"]

logging.config.dictConfig(LOGGING_CONFIG)


def select_by_date_modified(
    source: str | Path,
    year: int | None,
    month: int | None,
    day: int | None,
    *,
    suffixes: str = "nc",
    date: datetime,
) -> list[Path]:
    """Select files by the date on which they were last modified.

    Parameters
    ----------
    source : str or Path
    year : int
    month : int
    day : int
    suffixes : str
    date : datetime.date

    Returns
    -------
    list of Path
    """
    if date:
        date_selected = date
    else:
        date_selected = datetime(year, month, day)

    files = find_filepaths(source, file_suffixes=suffixes)

    selected_files = list()
    for file in files:
        if creation_date(file) == date_selected:
            logging.info(f"Selecting {file}.")
            selected_files.append(file)

    return selected_files
