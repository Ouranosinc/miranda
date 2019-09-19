import logging
from datetime import date
from datetime import datetime as dt
from pathlib import Path
from typing import List
from typing import Union

from miranda.utils import creation_date
from miranda.utils import find_filepaths

__all__ = ["select_by_date_modified"]


def select_by_date_modified(
    source: Union[Path, str],
    year: int,
    month: int,
    day: int,
    pattern: str = None,
    datetime: date = None,
) -> List:

    date_selected = date(year, month, day) or datetime

    files = find_filepaths(source, file_suffixes=pattern)

    selected_files = list()
    for file in files:
        if creation_date(file) == date_selected:
            logging.info(
                "{}: Selecting {}".format(dt.now().strftime("%Y-%m-%d %X"), file)
            )
            selected_files.append(file)

    return selected_files
