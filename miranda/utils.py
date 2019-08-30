import logging
import os
import platform
from contextlib import contextmanager
from datetime import date
from datetime import datetime as dt
from pathlib import Path
from types import GeneratorType
from typing import Iterable
from typing import List
from typing import Union

conversions = ["B", "k{}B", "M{}B", "G{}B", "T{}B", "P{}B", "E{}B", "Z{}B", "Y{}B"]


def size_formatter(i: int, binary: bool = True, precision: int = 2) -> str:
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


def file_size(
    file_path_or_bytes: Union[Path, str, int, List, GeneratorType],
    use_binary: bool = True,
    significant_digits: int = 2,
) -> str or None:
    """
    This function will return the size in bytes of a file or a list of files
    """

    if isinstance(file_path_or_bytes, int):
        return size_formatter(
            file_path_or_bytes, binary=use_binary, precision=significant_digits
        )
    elif isinstance(file_path_or_bytes, (list, GeneratorType)):
        sizes = [Path(f).stat().st_size for f in file_path_or_bytes]
        total = sum(sizes)
    elif Path(file_path_or_bytes).is_file():
        total = Path(file_path_or_bytes).stat().st_size
    else:
        return

    return size_formatter(total, binary=use_binary, precision=significant_digits)


def creation_date(path_to_file: Union[Path, str]) -> float or date:
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    if platform.system() == "Windows":
        return Path(path_to_file).stat().st_ctime
    else:
        stat = Path(path_to_file).stat()
        try:
            return date.fromtimestamp(stat.st_ctime)
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return date.fromtimestamp(stat.st_mtime)


def read_privileges(location: Union[Path, str]) -> bool:
    if os.access(location, os.R_OK):
        logging.info(
            "{}: {} Read OK!".format(dt.now().strftime("%Y-%m-%d %X"), location)
        )
    else:
        msg = "Ensure read privileges. Exiting."
        logging.error(msg)
        raise Exception(msg)
    return True


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


def find_files(
    source: Union[Path, str, GeneratorType, List[Union[Path, str]]],
    recursive: bool = True,
    file_suffixes: str = "*.nc",
    **_,
) -> (List, Path):

    if isinstance(source, (GeneratorType, List)):
        file_list = [Path(f) for f in source]
        source_path = os.path.commonpath(f for f in file_list)
    else:
        if recursive:
            file_list = [f for f in Path(source).rglob(file_suffixes)]
        elif not recursive:
            file_list = [f for f in Path(source).glob(file_suffixes)]
        else:
            raise ValueError("Recursive: {}".format(recursive))
        source_path = Path(source)
    return file_list, source_path


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
