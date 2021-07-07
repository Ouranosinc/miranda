import logging
import re
from collections import defaultdict
from logging import config
from pathlib import Path
from types import GeneratorType
from typing import List, Mapping, Union

from miranda.scripting import LOGGING_CONFIG
from miranda.storage import report_file_size
from miranda.utils import ingest

config.dictConfig(LOGGING_CONFIG)
Nested_List = List[List[Path]]
PathDict = Mapping[str, List[Path]]

GiB = int(pow(2, 30))

__all__ = [
    "group_by_deciphered_date",
    "group_by_length",
    "group_by_size",
    "group_by_subdirectories",
]


def group_by_length(files: Union[GeneratorType, List], size: int = 10) -> Nested_List:
    """
    This function groups files by an arbitrary number of file entries
    """
    logging.info("Creating groups of {} files".format(size))
    files = ingest(files)
    grouped_list = list()
    group = list()
    for i, f in enumerate(files):
        group.append(f)
        if (i + 1) % size == 0:
            grouped_list.append(group.copy())
            group.clear()
            continue
    if not group:
        pass
    else:
        grouped_list.append(group.copy())
    logging.info("Divided files into %s groups." % len(grouped_list))
    return grouped_list


def group_by_deciphered_date(files: Union[GeneratorType, List]) -> PathDict:
    """
    This function attempts to find a common date and groups files based on year and month
    """
    logging.info("Creating files from deciphered dates.")

    year_month_day = re.compile(
        r"(?P<year>[0-9]{4})-?(?P<month>[0-9]{2})-?(?P<day>[0-9]{2})?.*\.(?P<suffix>nc)$"
    )

    files = ingest(files)
    dates = defaultdict(lambda: list())
    total = 0
    for f in files:
        match = re.search(year_month_day, str(Path(f).name))
        if match.group("day"):
            key = "-".join([match.group("year"), match.group("month")])
            dates[key].append(Path(f))
            total += 1
        elif match.group("month"):
            key = match.group("year")
            dates[key].append(Path(f))
            total += 1
        else:
            continue

    if dates and total == len(files):
        logging.info(
            "All files have been grouped by date. {} groups created.".format(len(dates))
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
    files: Union[GeneratorType, List], size: int = 10 * GiB
) -> Nested_List:
    """
    This function will group files up until a desired size and save it as a grouping within a list
    """
    logging.info(
        "Creating groups of files based on size not exceeding {}.".format(
            report_file_size(size)
        )
    )

    files = ingest(files)
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
            continue
        elif total < size:
            continue

    if not group:
        logging.info("The final group is empty. Skipping this set...")
    else:
        grouped_list.append(group.copy())
    return grouped_list


def group_by_subdirectories(
    files: Union[GeneratorType, List], within: str or Path = None
) -> PathDict:
    """
    This function will group files based on the parent folder that they are located within.
    """
    if not within:
        within = Path.cwd()

    files = ingest(files)
    groups = defaultdict(list)
    for f in files:
        group_name = Path(f).relative_to(within).parent
        groups[group_name].append(f)

    logging.info(
        "File subdirectories found. Proceeding with {}.".format(
            str([str(key) for key in groups.keys()])
        )
    )
    return groups
