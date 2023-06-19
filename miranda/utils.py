"""Miscellaneous Helper Utilities module."""
from __future__ import annotations

import gzip
import logging.config
import os
import re
import sys
import tarfile
import tempfile
import warnings
import zipfile
from collections.abc import Iterable, Sequence
from contextlib import contextmanager
from io import StringIO
from pathlib import Path
from typing import TextIO

from .scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "HiddenPrints",
    "chunk_iterables",
    "generic_extract_archive",
    "list_paths_with_elements",
    "publish_release_notes",
    "single_item_list",
    "working_directory",
]

# For datetime validation
ISO_8601 = (
    r"^(-?(?:[1-9][0-9]*)?[0-9]{4})-(1[0-2]|0[1-9])-(3[01]|0[1-9]|[12][0-9])"
    r"T(2[0-3]|[01][0-9]):([0-5][0-9]):([0-5][0-9])(\.[0-9]+)?(Z|[+-](?:2[0-3]|[01][0-9]):[0-5][0-9])?$"
)


class HiddenPrints:
    """Special context manager for hiding print statements.

    Notes
    -----
    Solution from https://stackoverflow.com/a/45669280/7322852
    Credit to Alexander C (https://stackoverflow.com/users/2039471/alexander-c)
    CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)-
    """

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def chunk_iterables(iterable: Sequence, chunk_size: int) -> Iterable:
    """Generate lists of `chunk_size` elements from `iterable`.

    Notes
    -----
    Adapted from eidord (2012) https://stackoverflow.com/a/12797249/7322852 (https://creativecommons.org/licenses/by-sa/4.0/)
    """
    iterable = iter(iterable)
    while True:
        chunk = list()
        try:
            for _ in range(chunk_size):
                chunk.append(next(iterable))
            yield chunk
        except StopIteration:
            if chunk:
                yield chunk
            break


@contextmanager
def working_directory(directory: str | Path) -> None:
    """Change the working directory within a context object.

    This function momentarily changes the working directory within the context and reverts to the file working directory
    when the code block it is acting upon exits

    Parameters
    ----------
    directory : str or pathlib.Path

    Returns
    -------
    None
    """
    owd = os.getcwd()

    if (2, 7) < sys.version_info < (3, 6):
        directory = str(directory)

    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)


def single_item_list(iterable: Iterable) -> bool:
    """Ascertain whether a list has exactly one entry.

    See: https://stackoverflow.com/a/16801605/7322852

    Parameters
    ----------
    iterable : Iterable

    Returns
    -------
    bool

    """
    iterator = iter(iterable)
    has_true = any(iterator)  # consume from "i" until first true or it's exhausted
    has_another_true = any(
        iterator
    )  # carry on consuming until another true value / exhausted
    return has_true and not has_another_true  # True if exactly one true found


def generic_extract_archive(
    resources: str | Path | list[bytes | str | Path],
    output_dir: str | Path | None = None,
) -> list[Path]:
    """Extract archives (tar/zip) to a working directory.

    Parameters
    ----------
    resources : str or Path or list of bytes or str or Path
        list of archive files (if netCDF files are in list, they are passed and returned as well in the return).
    output_dir : str or Path, optional
        string or Path to a working location (default: temporary folder).

    Returns
    -------
    list
        List of original or of extracted files
    """
    archive_types = [".gz", ".tar", ".zip", ".7z"]
    output_dir = output_dir or tempfile.gettempdir()

    if not isinstance(resources, list):
        resources = [resources]

    files = list()

    for arch in resources:
        if any(ext in str(arch).lower() for ext in archive_types):
            try:
                logging.debug("archive=%s", arch)
                file = Path(arch)

                if file.suffix == ".nc":
                    files.append(Path(output_dir.join(arch)))
                elif file.suffix == ".tar":
                    with tarfile.open(arch, mode="r") as tar:
                        tar.extractall(path=output_dir)
                        files.extend(
                            [Path(output_dir).joinpath(f) for f in tar.getnames()]
                        )
                elif file.suffix == ".zip":
                    with zipfile.ZipFile(arch, mode="r") as zf:
                        zf.extractall(path=output_dir)
                        files.extend(
                            [Path(output_dir).joinpath(f) for f in zf.namelist()]
                        )
                elif file.suffix == ".gz":
                    logging.warning(
                        "GZIP file found. Can only extract one expected file."
                    )
                    with gzip.open(arch, "rb") as gf, open(
                        Path(output_dir).joinpath(arch.stem), "w"
                    ) as f_out:
                        f_out.write(gf.read().decode("utf-8"))
                        files.append(Path(output_dir).joinpath(arch.stem))
                elif file.suffix == ".7z":
                    msg = "7z file extraction is not supported at this time."
                    logging.warning(msg)
                    warnings.warn(msg, UserWarning)
                else:
                    logging.debug('File extension "%s" unknown' % file)
            except Exception as e:
                logging.error(f"Failed to extract sub archive {arch}: {e}")
        else:
            logging.warning("No archives found. Continuing...")
            return resources

    return files


########################################################################################


def list_paths_with_elements(
    base_paths: str | list[str], elements: list[str]
) -> list[dict]:
    """List a given path structure.

    Parameters
    ----------
    base_paths : list of str
        List of paths from which to start the search.
    elements : list of str
        Ordered list of the expected elements.

    Returns
    -------
    list of dict
        The keys are 'path' and each of the members of the given elements, the path is the absolute path.

    Notes
    -----
    Suppose you have the following structure: /base_path/{color}/{shape}
    The resulting list would look like::

         [{'path':/base_path/red/square, 'color':'red', 'shape':'square'},
         {'path':/base_path/red/circle, 'color':'red', 'shape':'circle'},
         {'path':/base_path/blue/triangle, 'color':'blue', 'shape':'triangle'},
         ...]

    Obviously, 'path' should not be in the input list of elements.
    """
    # Make sure the base_paths input is a list of absolute path
    paths = list()
    if not hasattr(base_paths, "__iter__"):
        paths.append(base_paths)
    paths = map(os.path.abspath, base_paths)
    # If elements list is empty, return empty list (end of recursion).
    if not elements:
        return list()

    paths_elements = list()
    for base_path in paths:
        try:
            path_content = [f for f in Path(base_path).iterdir()]
        except NotADirectoryError:
            continue
        path_content.sort()
        next_base_paths = []
        for path_item in path_content:
            next_base_paths.append(base_path.joinpath(path_item))
        next_pe = list_paths_with_elements(next_base_paths, elements[1:])
        if next_pe:
            for i, one_pe in enumerate(next_pe):
                relative_path = next_pe[i]["path"].replace(base_path, "", 1)
                new_element = relative_path.split("/")[1]
                next_pe[i][elements[0]] = new_element
            paths_elements.extend(next_pe)
        elif len(elements) == 1:
            for my_path, my_item in zip(next_base_paths, path_content):
                paths_elements.append({"path": my_path, elements[0]: my_item})
    return paths_elements


def publish_release_notes(
    style: str = "md", file: os.PathLike | StringIO | TextIO | None = None
) -> str | None:
    """Format release history in Markdown or ReStructuredText.

    Parameters
    ----------
    style : {"rst", "md"}
        Use ReStructuredText formatting or Markdown. Default: Markdown.
    file : {os.PathLike, StringIO, TextIO}, optional
        If provided, prints to the given file-like object.
        Otherwise, returns a string.

    Returns
    -------
    str, optional

    Notes
    -----
    This function is solely for development purposes.
    """
    history_file = Path(__file__).parent.parent.joinpath("HISTORY.rst")

    if not history_file.exists():
        raise FileNotFoundError("History file not found in miranda file tree.")

    with open(history_file) as hf:
        history = hf.read()

    if style == "rst":
        hyperlink_replacements = {
            r":issue:`([0-9]+)`": r"`GH/\1 <https://github.com/Ouranosinc/miranda/issues/\1>`_",
            r":pull:`([0-9]+)`": r"`PR/\1 <https://github.com/Ouranosinc/miranda/pull/\>`_",
            r":user:`([a-zA-Z0-9_.-]+)`": r"`@\1 <https://github.com/\1>`_",
        }
    elif style == "md":
        hyperlink_replacements = {
            r":issue:`([0-9]+)`": r"[GH/\1](https://github.com/Ouranosinc/miranda/issues/\1)",
            r":pull:`([0-9]+)`": r"[PR/\1](https://github.com/Ouranosinc/miranda/pull/\1)",
            r":user:`([a-zA-Z0-9_.-]+)`": r"[@\1](https://github.com/\1)",
        }
    else:
        raise NotImplementedError()

    for search, replacement in hyperlink_replacements.items():
        history = re.sub(search, replacement, history)

    if style == "md":
        history = history.replace(".. :changelog:\n\n", "")
        history = history.replace("=======\nHistory\n=======", "# History")

        titles = {r"\n(.*?)\n([\-]{1,})": "-", r"\n(.*?)\n([\^]{1,})": "^"}
        for title_expression, level in titles.items():
            found = re.findall(title_expression, history)
            for grouping in found:
                fixed_grouping = (
                    str(grouping[0]).replace("(", r"\(").replace(")", r"\)")
                )
                search = rf"({fixed_grouping})\n([\{level}]{'{' + str(len(grouping[1])) + '}'})"
                replacement = f"{'##' if level=='-' else '###'} {grouping[0]}"
                history = re.sub(search, replacement, history)

        link_expressions = r"[\`]{1}([\w\s]+)\s<(.+)>`\_"
        found = re.findall(link_expressions, history)
        for grouping in found:
            search = rf"`{grouping[0]} <.+>`\_"
            replacement = f"[{str(grouping[0]).strip()}]({grouping[1]})"
            history = re.sub(search, replacement, history)

    if not file:
        return history
    if isinstance(file, (Path, os.PathLike)):
        file = Path(file).open("w")
    print(history, file=file)


def read_privileges(location: str | Path, strict: bool = False) -> bool:
    """Determine whether a user has read privileges to a specific file.

    Parameters
    ----------
    location: str or Path
    strict: bool

    Returns
    -------
    bool
        Whether the current user shell has read privileges
    """
    msg = ""
    try:
        if Path(location).exists():
            if os.access(location, os.R_OK):
                msg = f"{location} is read OK!"
                logging.info(msg)
                return True
            msg = f"Ensure read privileges for `{location}`."
        else:
            msg = f"`{location}` is an invalid path."
        raise OSError()

    except OSError:
        logging.exception(msg)
        if strict:
            raise
        return False
