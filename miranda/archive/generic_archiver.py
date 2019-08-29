#!/bin/env python3
import logging
import os
import re
from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime as dt
from math import pow
from pathlib import Path
from types import GeneratorType
from typing import Iterable
from typing import List
from typing import Mapping

Nested_List = List[List[Path]]
PathDict = Mapping[str, List[Path]]

GiB = int(pow(2, 30))


@contextmanager
def working_directory(directory):
    """
    This function momentarily changes the working directory within the
     context and reverts to the file working directory when the code block
     it is acting upon exits
    """
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)
    return


def file_size(
    file_path_or_bytes: str or Path or int,
    use_binary: bool = True,
    significant_digits: int = 2,
) -> str or None:
    """
    This function will return the size in bytes of a file or a list of files
    """

    conversions = ["B", "k{}B", "M{}B", "G{}B", "T{}B", "P{}B", "E{}B", "Z{}B", "Y{}B"]

    def _size_formatter(i: int, binary: bool = True, precision: int = 2) -> str:
        """
        This function will format byte size into an appropriate nomenclature
        """
        import math

        base = 1024 if binary else 1000
        if i == 0:
            return "0 B"
        multiple = math.trunc(math.log2(i) / math.log2(base))
        value = i / math.pow(base, multiple)
        suffix = conversions[multiple].format("i" if binary else "")
        return "{value:.{precision}f} {suffix}".format(**locals())

    if isinstance(file_path_or_bytes, int):
        return _size_formatter(
            file_path_or_bytes, binary=use_binary, precision=significant_digits
        )
    elif isinstance(file_path_or_bytes, (list, GeneratorType)):
        sizes = [Path(f).stat().st_size for f in file_path_or_bytes]
        total = sum(sizes)
    elif Path(file_path_or_bytes).is_file():
        total = Path(file_path_or_bytes).stat().st_size
    else:
        return

    return _size_formatter(total, binary=use_binary, precision=significant_digits)


def single_item_list(iterable: Iterable) -> bool:
    """
    See: https://stackoverflow.com/a/16801605/7322852
    """
    iterator = iter(iterable)
    has_true = any(iterator)  # consume from "i" until first true or it's exhausted
    has_another_true = any(
        iterator
    )  # carry on consuming until another true value / exhausted
    return has_true and not has_another_true  # True if exactly one true found


def group_by_length(files: list or GeneratorType, size: int = 10) -> Nested_List:
    """
    This function groups files by an arbitrary number of file entries
    """

    grouped_list = list()
    group = list()
    if isinstance(files, GeneratorType):
        files = [f for f in files]
    files.sort()

    logging.info(
        "{}: Creating groups of {} files".format(dt.now().strftime("%Y-%m-%d %X"), size)
    )

    for i, f in enumerate(files):
        group.append(f)
        if (i + 1) % size == 0:
            grouped_list.append(group.copy())
            group.clear()
            continue

    if not group:
        logging.info(
            "{}: The final group is empty. Skipping...".format(
                dt.now().strftime("%Y-%m-%d %X")
            )
        )
    else:
        grouped_list.append(group.copy())
    return grouped_list


def group_by_deciphered_date(files: list or GeneratorType) -> PathDict:
    """
    This function attempts to find a common date and groups files based on year and month
    """
    if isinstance(files, GeneratorType):
        files = [Path(f) for f in files]
    files.sort()

    logging.info(
        "{}: Creating files from deciphered dates.".format(
            dt.now().strftime("%Y-%m-%d %X")
        )
    )

    year_month_day = re.compile(
        r"(?P<year>[0-9]{4})-?(?P<month>[0-9]{2})-?(?P<day>[0-9]{2})?.*\.(?P<suffix>nc)$"
    )

    dates = defaultdict(lambda: list())
    total = 0
    for f in files:
        match = re.search(year_month_day, str(f.name))
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

    now = dt.now()
    if dates and total == len(files):
        logging.info(
            "{}: All files have been grouped by date.".format(
                now.strftime("%Y-%m-%d %X")
            )
        )
        return dates

    elif dates and total != len(files):
        logging.info(
            "{}: Not all files were successfully grouped by date. Grouping aborted.".format(
                now.strftime("%Y-%m-%d %X")
            )
        )
    else:
        logging.info(
            "{}: No matches for dates found. Grouping aborted.".format(
                now.strftime("%Y-%m-%d %X")
            )
        )
    return dict(data=files)


def group_by_size(files: list or GeneratorType, size: int = 10 * GiB) -> Nested_List:
    """
    This function will group files up until a desired size and save it as a grouping within a list
    """
    grouped_list = list()
    group = list()
    if isinstance(files, GeneratorType):
        files = [f for f in files]
    files.sort()

    logging.info(
        "{}: Creating groups of files based on size not exceeding {}".format(
            dt.now().strftime("%Y-%m-%d %X"), file_size(size)
        )
    )

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
        logging.info(
            "{}: The final group is empty. Skipping...".format(
                dt.now().strftime("%Y-%m-%d %X")
            )
        )
    else:
        grouped_list.append(group.copy())
    return grouped_list


def group_by_subdirectories(
    files: list or GeneratorType, within: str or Path = None
) -> PathDict:
    """
    This function will group files based on the parent folder that they are located within.
    """
    groupings = defaultdict(list)
    if isinstance(files, GeneratorType):
        files = [f for f in files]
    files.sort()

    if not within:
        within = Path.cwd()

    for f in files:
        group_name = Path(f).relative_to(within).parent
        groupings[group_name].append(f)

    logging.info(
        "{}: File subdirectories found. Proceeding with {}.".format(
            dt.now().strftime("%Y-%m-%d %X"),
            str([str(key) for key in groupings.keys()]),
        )
    )
    return groupings


if __name__ == "__main__":
    logging.basicConfig(
        filename="{}_{}.log".format(
            dt.strftime(dt.now(), "%Y%m%d"), Path(__name__).stem
        ),
        level=logging.INFO,
    )
