import platform
from datetime import date
from pathlib import Path
from types import GeneratorType

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
    file_path_or_bytes: str or Path or int,
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


def creation_date(path_to_file: str or Path) -> float or date:
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
