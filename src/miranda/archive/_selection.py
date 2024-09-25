"""Select files by the date on which they were last modified."""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

from miranda.io import find_filepaths
from miranda.io.utils import creation_date
from miranda.scripting import LOGGING_CONFIG

__all__ = ["select_by_date_modified"]

logging.config.dictConfig(LOGGING_CONFIG)  # noqa


def select_by_date_modified(
    source: str | Path,
    year: int | None,
    month: int | None,
    day: int | None,
    *,
    suffixes: str = "nc",
    date: datetime.date,
) -> list[Path]:
    """
    Select files by the date on which they were last modified.

    Parameters
    ----------
    source : str or Path
        The directory to search for files.
    year : int
        The year of the date to select.
    month : int
        The month of the date to select
    day : int
        The day of the date to select.
    suffixes : str
        The file suffixes to search.
    date : date
        The date to select.

    Returns
    -------
    list of Path
        The selected files.
    """
    if date:
        date_selected = date
    else:
        date_selected = datetime(year, month, day)

    files = find_filepaths(source, file_suffixes=suffixes)

    selected_files = list()
    for file in files:
        if creation_date(file) == date_selected:
            msg = f"Selecting {file}."
            logging.info(msg)
            selected_files.append(file)

    return selected_files
