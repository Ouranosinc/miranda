"""Miscellaneous Helper Utilities module."""

from __future__ import annotations
import gzip
import logging
import os
import sys
import tarfile
import tempfile
import warnings
import zipfile
from collections.abc import Iterable, Sequence
from contextlib import contextmanager
from pathlib import Path


__all__ = [
    "HiddenPrints",
    "chunk_iterables",
    "generic_extract_archive",
    "list_paths_with_elements",
    "single_item_list",
    "working_directory",
]

from types import GeneratorType


# For datetime validation
ISO_8601 = (
    r"^(-?(?:[1-9][0-9]*)?[0-9]{4})-(1[0-2]|0[1-9])-(3[01]|0[1-9]|[12][0-9])"
    r"T(2[0-3]|[01][0-9]):([0-5][0-9]):([0-5][0-9])(\.[0-9]+)?(Z|[+-](?:2[0-3]|[01][0-9]):[0-5][0-9])?$"
)


class HiddenPrints:
    """
    Special context manager for hiding print statements.

    Notes
    -----
    Solution from https://stackoverflow.com/a/45669280/7322852
    Credit to Alexander C (https://stackoverflow.com/users/2039471/alexander-c)
    CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)-
    """

    def __enter__(self):  # noqa: D105
        self._original_stdout = sys.stdout
        sys.stdout = Path(os.devnull).open("w")

    def __exit__(self, exc_type, exc_val, exc_tb):  # noqa: D105
        sys.stdout.close()
        sys.stdout = self._original_stdout


def chunk_iterables(iterable: Sequence, chunk_size: int) -> Iterable:
    """
    Generate lists of `chunk_size` elements from `iterable`.

    Parameters
    ----------
    iterable : Sequence
        The iterable to chunk.
    chunk_size : int
        The size of the chunks.

    Yields
    ------
    Iterable
        The chunked iterable.

    Notes
    -----
    Adapted from eidord (2012) https://stackoverflow.com/a/12797249/7322852 (https://creativecommons.org/licenses/by-sa/4.0/)
    """
    iterable = iter(iterable)
    while True:
        chunk = list()
        try:
            for _ in range(chunk_size):
                chunk.append(next(iterable))  # noqa: PERF401
            yield chunk
        except StopIteration:
            if chunk:
                yield chunk
            break


# FIXME: The following function could probably be replaced or at least placed closer to its usages.
@contextmanager
def working_directory(directory: str | Path) -> None:
    """
    Change the working directory within a context object.

    This function momentarily changes the working directory within the context and reverts to the file working directory
    when the code block it is acting upon exits

    Parameters
    ----------
    directory : str or pathlib.Path
        The directory to temporarily change to.
    """
    owd = os.getcwd()  # noqa: PTH109

    if (2, 7) < sys.version_info < (3, 6):
        directory = str(directory)

    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)


# FIXME: The following function could probably be replaced or at least placed closer to its usages.
def group_by_length(
    files: GeneratorType | list[str | Path],
    size: int = 10,
    sort: bool = False,
) -> list[list[Path]]:
    """
    Group files by an arbitrary number of file entries.

    Parameters
    ----------
    files : GeneratorType or list of str or pathlib.Path
        The files to be grouped.
    size : int
        The number of files to be grouped together.
    sort : bool
        Sort the files before grouping.

    Returns
    -------
    list[list[pathlib.Path]]
        Grouped files.
    """
    msg = f"Creating groups of {size} files"
    logging.info(msg)
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
    msg = f"Divided files into {len(grouped_list)} groups."
    logging.info(msg)
    return grouped_list


# FIXME: The following function could probably be replaced or at least placed closer to its usages.
def single_item_list(iterable: Iterable) -> bool:
    """
    Ascertain whether a list has exactly one entry.

    See: https://stackoverflow.com/a/16801605/7322852

    Parameters
    ----------
    iterable : Iterable
        The list to check.

    Returns
    -------
    bool
        Whether the list is a single item.
    """
    iterator = iter(iterable)
    has_true = any(iterator)  # consume from "i" until first true or it's exhausted
    has_another_true = any(iterator)  # carry on consuming until another true value / exhausted
    return has_true and not has_another_true  # True if exactly one true found


def generic_extract_archive(
    resources: str | Path | list[bytes | str | Path],
    output_dir: str | Path | None = None,
) -> list[Path]:
    """
    Extract archives (tar/zip) to a working directory.

    Parameters
    ----------
    resources : str or Path or list of bytes or str or Path
        List of archive files (if netCDF files are in list, they are passed and returned as well in the return).
    output_dir : str or Path, optional
        String or Path to a working location (default: temporary folder).

    Returns
    -------
    list
        The list of original or of extracted files.
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
                        safe_extract(tar, path=output_dir)
                        files.extend([Path(output_dir).joinpath(f) for f in tar.getnames()])
                elif file.suffix == ".zip":
                    with zipfile.ZipFile(arch, mode="r") as zf:
                        safe_extract(zf, path=output_dir)
                        files.extend([Path(output_dir).joinpath(f) for f in zf.namelist()])
                elif file.suffix == ".gz":
                    logging.warning("GZIP file found. Can only extract one expected file.")
                    with (
                        gzip.open(arch, "rb") as gf,
                        Path(output_dir).joinpath(arch.stem).open("w") as f_out,
                    ):
                        f_out.write(gf.read().decode("utf-8"))
                        files.append(Path(output_dir).joinpath(arch.stem))
                elif file.suffix == ".7z":
                    msg = "7z file extraction is not supported at this time."
                    logging.warning(msg)
                    warnings.warn(msg, UserWarning, stacklevel=2)
                else:
                    msg = f'File extension "{file}") unknown'
                    logging.debug(msg)
            except OSError as e:
                msg = f"Failed to extract sub archive {arch}: {e}"
                logging.error(msg)
        else:
            logging.warning("No archives found. Continuing...")
            return resources

    return files


########################################################################################


def list_paths_with_elements(base_paths: str | list[str] | os.PathLike[str], elements: list[str]) -> list[dict]:
    """
    List a given path structure.

    Parameters
    ----------
    base_paths : str or list of str or os.PathLike
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
    paths = []
    if not hasattr(base_paths, "__iter__"):
        paths.append(base_paths)
    paths = map(os.path.abspath, base_paths)
    # If elements list is empty, return empty list (end of recursion).
    if not elements:
        return []

    paths_elements = []
    for base_path in paths:
        try:
            path_content = [f for f in Path(base_path).iterdir()]
        except NotADirectoryError:
            msg = "Not a directory. Skipping..."
            logging.debug(msg)
            continue
        path_content.sort()
        next_base_paths = [Path(base_path).joinpath(path_item) for path_item in path_content]
        next_pe = list_paths_with_elements(next_base_paths, elements[1:])
        if next_pe:
            for i in range(len(next_pe)):
                relative_path = next_pe[i]["path"].replace(base_path, "", 1)
                new_element = relative_path.split("/")[1]
                next_pe[i][elements[0]] = new_element
            paths_elements.extend(next_pe)
        elif len(elements) == 1:
            for my_path, my_item in zip(next_base_paths, path_content, strict=False):
                paths_elements.append({"path": my_path, elements[0]: my_item})
    return paths_elements


def read_privileges(location: str | Path, strict: bool = False) -> bool:
    """
    Determine whether a user has read privileges to a specific file.

    Parameters
    ----------
    location : str or Path
        The location to be assessed.
    strict : bool
        Whether to raise an exception if the user does not have read privileges. Default: False.

    Returns
    -------
    bool
        Whether the current user shell has read privileges.
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


def _is_within_directory(directory: str | os.PathLike, target: str | os.PathLike) -> bool:
    """
    Check if a target path is within a directory.

    Parameters
    ----------
    directory : str or os.PathLike
        The directory to check.
    target : str or os.PathLike
        The target path to check.

    Returns
    -------
    bool
        Whether the target path is within the directory.

    Notes
    -----
    Function addressing exploit CVE-2007-4559 for both tar and zip files.
    """
    abs_directory = Path(directory).resolve()
    abs_target = Path(target).resolve()

    prefix = os.path.commonprefix([abs_directory, abs_target])
    return prefix == abs_directory


def safe_extract(
    archive: tarfile.TarFile | zipfile.ZipFile,
    path: str = ".",
    members: list[str] | None = None,
    *,
    numeric_owner: bool = False,
) -> None:
    """
    Extract all members from the archive to the current working directory or directory path.

    Parameters
    ----------
    archive : TarFile or ZipFile
        The archive to extract.
    path : str, optional
        The path to extract the archive to.
    members : list of str, optional
        The members to extract.
    numeric_owner : bool
        Whether to extract the archive with numeric owner. Default: False.

    Notes
    -----
    Function addressing exploit CVE-2007-4559 for both tar and zip files.
    """
    if isinstance(archive, tarfile.TarFile):
        for member in archive.getmembers():
            member_path = Path(path).joinpath(member.name)
            if not _is_within_directory(path, member_path):
                raise Exception("Attempted Path Traversal in Tar File")
        archive.extractall(  # noqa: S202
            path, members=members, numeric_owner=numeric_owner
        )

    elif isinstance(archive, zipfile.ZipFile):
        for member in archive.namelist():
            member_path = Path(path).joinpath(member)
            if not _is_within_directory(path, member_path):
                raise Exception("Attempted Path Traversal in Zip File")
        archive.extractall(path, members=members)  # noqa: S202

    else:
        raise TypeError("Archive must be a TarFile or ZipFile object.")
