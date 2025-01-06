from __future__ import annotations

import logging
import re
from logging.config import dictConfig
from pathlib import Path
from types import GeneratorType
from typing import Dict, List

from miranda.scripting import LOGGING_CONFIG
from miranda.storage import report_file_size

dictConfig(LOGGING_CONFIG)
Nested_List = List[List[Path]]
PathDict = Dict[str, List[Path]]


GiB = int(pow(2, 30))

__all__ = [
    "group_by_deciphered_date",
    "group_by_length",
    "group_by_size",
    "group_by_subdirectories",
]


def group_by_length(
    files: GeneratorType | list[str | Path],
    size: int = 10,
    sort: bool = False,
) -> list[list[Path]]:
    """Group files by an arbitrary number of file entries.

    Parameters
    ----------
    files: GeneratorType or list of str or pathlib.Path
    size: int
    sort: bool

    Returns
    -------
    list[list[pathlib.Path]]
    """
    logging.info(f"Creating groups of {size} files")
    if sort:
        files = [Path(f) for f in files]
        files.sort()
    grouped_list = list()
    group = list()
    for i, f in enumerate(files):
        group.append(Path(f))
        if (i + 1) % size == 0:
            grouped_list.append(group.copy())
            group.clear()
            continue
    if not group:
        pass
    else:
        grouped_list.append(group.copy())
    logging.info(f"Divided files into {len(grouped_list)} groups.")
    return grouped_list


def group_by_deciphered_date(
    files: GeneratorType | list[str | Path],
) -> dict[str, list[Path]]:
    """Find a common date and groups files based on year and month.

    Parameters
    ----------
    files: GeneratorType or list of str or pathlib.Path

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    logging.warning("This function doesn't work well with multi-thread processing!")
    logging.info("Creating files from deciphered dates.")

    year_month_day = re.compile(
        r"(?P<year>\d{4})-?(?P<month>\d{2})-?(?P<day>\d{2})?.*\.(?P<suffix>nc|zarr)$"
    )

    files = [Path(f) for f in files]
    files.sort()
    dates = dict()
    total = 0
    for f in files:
        match = re.search(year_month_day, str(Path(f).name))
        if match.group("day"):
            key = "-".join([match.group("year"), match.group("month")])
            dates.setdefault(key, list()).append(Path(f))
            total += 1
        elif match.group("month"):
            key = match.group("year")
            dates.setdefault(key, list()).append(Path(f))
            total += 1
        else:
            continue

    if dates and total == len(files):
        logging.info(
            f"All files have been grouped by date. {len(dates)} groups created."
        )
        return dict(dates)

    if dates and total != len(files):
        logging.info(
            "Not all files were successfully grouped by date. Grouping aborted."
        )
    else:
        logging.info("No matches for dates found. Grouping aborted.")
    return dict(data=files)


def group_by_size(
    files: GeneratorType | list[str | Path], size: int = 10 * GiB
) -> list[list[Path]]:
    """Group files up until a desired size and save it as a grouping within a list.

    Parameters
    ----------
    files : GeneratorType or list of str or pathlib.Path
    size : int

    Returns
    -------
    list[list[pathlib.Path]]
    """
    logging.info(
        f"Creating groups of files based on size not exceeding: {report_file_size(size)}."
    )

    files = [Path(f) for f in files]
    files.sort()
    grouped_list = list()
    group = list()
    total = 0
    for f in files:
        total += Path.stat(f).st_size
        group.append(f)
        if total > size:
            grouped_list.append(group.copy())
            group.clear()
            total = 0

    if not group:
        logging.info("The final group is empty. Skipping this set...")
    else:
        grouped_list.append(group.copy())
    return grouped_list


def group_by_subdirectories(
    files: GeneratorType | list[str | Path], within: str | Path = None
) -> dict[str, list[Path]]:
    """Group files based on the parent folder that they are located within.

    Parameters
    ----------
    files : GeneratorType or list of str or pathlib.Path
    within : str or pathlib.Path

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    if not within:
        within = Path.cwd()

    files = [Path(f) for f in files]
    files.sort()
    groups = dict()
    for f in files:
        group_name = Path(f).relative_to(within).parent
        groups.setdefault(group_name, list()).append(f)

    logging.info(
        f"File subdirectories found. Proceeding with: `{', '.join([str(key) for key in groups.keys()])}`."
    )
    return groups
