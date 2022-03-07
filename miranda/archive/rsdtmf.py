"""
=================
rstdmf Operations
=================

Classes:

 * rstdmfError - the exception raised on failure.

Objects:

 * DiskSpaceEvent - threading event for when there is not enough space on disk.

Functions:

 * :func:`rstdmf_rename` - single rstdmf call with file renaming at the end.
 * :func:`local_storage_for_rstdmf` - storage object for rstdmf_divisions.
 * :func:`rstdmf_divisions` - rstdmf calls based on file numbers and size.

Notes
-----
The restored files are always renamed to their full path with '/' replaced
by '.' (leading '/' removed) to allow for restoring files with the same name
in different locations.
e.g. /server/path/sample.tar becomes server.path.sample.tar
The python trick to convert is::

   >>> path_string = Path(dmf1_file).absolute() # to get full path if necessary.
   >>> str(path_string).replace('/','.').lstrip('.')

"""
import logging
import os
import subprocess
import threading
import time
from logging import config
from pathlib import Path
from typing import List

from miranda.archive.ops import transfer_file
from miranda.scripting import LOGGING_CONFIG
from miranda.storage import (
    DiskSpaceError,
    FileMeta,
    StorageState,
    size_division,
    size_evaluation,
)
from miranda.units import GiB
from miranda.utils import find_filepaths, yesno_prompt

DiskSpaceEvent = threading.Event()

config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "local_storage_for_rstdmf",
    "rstdmf_divisions",
    "rstdmf_rename",
    "RstdmfError",
]


def rstdmf_rename(file_list: List[str], restore_path: str):
    """Single rstdmf call with file renaming at the end.

    Parameters
    ----------
    file_list : List[str]
      list of files to restore.
    restore_path : str

    Notes
    -----
    The restored files are renamed to their full path with '/' replaced by '.'.
        e.g. /dmf1/scenario/sample.tar becomes dmf1.scenario.sample.tar

    """

    # Make sure we have the full path of each argument
    if not hasattr(file_list, "__iter__"):
        file_list = [file_list]
    file_list = map(Path().absolute(), file_list)
    restore_path = Path(restore_path).absolute()
    # Create output directory, if necessary
    try:
        Path(restore_path).mkdir(parents=True, exist_ok=True)
    except OSError:
        raise RstdmfError(f"Cannot create restore path {restore_path}.")
    # Create string with list of files
    file_list_string = "' '".join(file_list)
    file_list_string = "'" + file_list_string + "'"
    file_list_string = file_list_string.replace("(", r"\(")
    file_list_string = file_list_string.replace(")", r"\)")
    # rstdmf call
    flag_system = subprocess.call(
        ["rstdmf", "-m", "-b", file_list_string, restore_path]
    )
    if flag_system:
        raise RstdmfError("rstdmf call failed.")
    # rename files
    for file_to_restore in file_list:
        file_name = file_to_restore.replace("/", ".").lstrip(".")
        restored_file = Path(restore_path).joinpath(Path(file_to_restore).name)
        try:
            transfer_file(restored_file, Path(restore_path).joinpath(file_name))
        except OSError:
            raise RstdmfError(f"Could not move {restored_file}.")


def local_storage_for_rstdmf(
    files_for_rstdmf: List[str], restore_path: str, max_size_on_disk: int
):
    """Create a local storage object for rstdmf_divisions function.

    Parameters
    ----------
    files_for_rstdmf : List[str]
        files that will be restored.
    restore_path : str
    max_size_on_disk : int
        maximum amount of space to be used on disk.

    Returns
    -------
    out : StorageState object

    """

    # Get files already restored on disk and their size
    files_on_disk = []
    for file_on_disk in files_for_rstdmf:
        renamed_file = file_on_disk.replace("/", ".").lstrip(".")
        renamed_path = Path(restore_path).joinpath(renamed_file)
        if os.path.isfile(renamed_path):
            files_on_disk.append(renamed_path)
        elif os.path.isdir(renamed_path):
            files_within = find_filepaths(renamed_path)
            files_on_disk.append(
                [FileMeta(renamed_path), size_evaluation(files_within)]
            )
    size_on_disk = size_evaluation(files_on_disk)
    # create a local storage object
    return StorageState(
        restore_path, max_size_on_disk, size_on_disk, max_size_on_disk - size_on_disk
    )


#
def rstdmf_divisions(
    paths_to_restore,
    restore_path,
    disk_space_margin: int = 100 * GiB,
    disk_refresh_time: int = 60,
    rstdmf_size_limit: int = 10 * GiB,
    rstdmf_file_limit: int = 24,
    max_size_on_disk: int = 0,
    verbose: bool = False,
    preserve_order: bool = True,
):
    """Securely make call(s) to rstdmf.

    Parameters
    ----------
    paths_to_restore : list of strings
        list of files to restore (path or FileMeta objects), giving
        directories also works but they are considered like a single file
        (i.e. be careful of their size!).
    restore_path : str
        Path where files will be restored.
    disk_space_margin : int
        minimum space required on disk to proceed, set to 0 for no minimum
        requirement (default: 100 GiB).
    disk_refresh_time : int
        time to wait to lookup disk space again when full (default: 60 s).
    rstdmf_size_limit : int
        size limit of a single rstdmf call, set to 0 for no size limit
        (default: 10 GiB).
    rstdmf_file_limit : int
        limit in number of files of a single rstdmf call, set to 0 for no
        limit (default: 24).
    max_size_on_disk : int
        maximum size of restored files on disk, set to 0 for no maximum size
        (default: 0).
    verbose : bool
        verbose mode (default: off).
    preserve_order : bool
        flag to force files to be restored in the order they are given
        (default: True).

    Notes
    -----
    The restored files are renamed to their full path with '/' replaced by '.':
    for example, ``/server/path/sample.tar`` becomes ``server.path.sample.tar``
    This allows getting files with the same name in different locations.
    If this renamed file is already in the target location, it is not restored.

    """

    # Make a list of files that should be restored
    files_for_rstdmf = []
    for file_to_restore in paths_to_restore:
        # If a directory is given, create its associated FileMeta object
        if os.path.isdir(file_to_restore):
            files_within = find_filepaths(file_to_restore)
            file_to_restore = FileMeta(file_to_restore, size_evaluation(files_within))
        # If just a file name is given, convert it to a FileMeta object
        elif not isinstance(file_to_restore, FileMeta):
            file_to_restore = FileMeta(file_to_restore)
        # Skip files that do not exist or are already in the restoration path
        file_name = file_to_restore.path.replace("/", ".").lstrip(".")
        if (
            Path(restore_path).joinpath(file_name).exists()
            or not Path(file_to_restore).exists()
        ):
            continue
        # Make sure no files are repeated
        for file_for_rstdmf in files_for_rstdmf:
            if file_for_rstdmf.path == file_to_restore.path:
                break
        else:
            files_for_rstdmf.append(file_to_restore)
    if not files_for_rstdmf:
        logging.warning("All files are already on disk, or nonexistant.")
        return
    # Divide files to satisfy rstdmf constraints
    if (max_size_on_disk != 0) and (max_size_on_disk < rstdmf_size_limit):
        division_size_limit = max_size_on_disk
    else:
        division_size_limit = rstdmf_size_limit
        divisions = size_division(
            files_for_rstdmf,
            division_size_limit,
            rstdmf_file_limit,
            check_name_repetition=True,
            preserve_order=preserve_order,
        )
    try:
        storage = StorageState(restore_path)
    except DiskSpaceError:
        storage = StorageState(restore_path, 0, 0, 0)
    if verbose:
        size_of_files = size_evaluation(files_for_rstdmf)
        logging.warning(
            "Disk space: {} GB, Restoration size: {} GiB, ".format(
                str(float(storage.free_space) / GiB), str(size_of_files / GiB)
            )
            + f"Subdivisions: {str(len(divisions))}"
        )
        if not yesno_prompt("Proceed with restoration?"):
            return

    for division in divisions:
        size_of_files = size_evaluation(division)
        list_of_path = [filemeta.path for filemeta in division]
        # If we can't get disk space, disregard it in the process
        try:
            storage = StorageState(restore_path)
        except DiskSpaceError:
            storage = StorageState(restore_path, 0, 0, 0)
            disk_space_margin = 0
        local_storage = local_storage_for_rstdmf(
            files_for_rstdmf, restore_path, max_size_on_disk
        )
        if (
            (disk_space_margin != 0)
            and (storage.free_space - disk_space_margin < size_of_files)
        ) or ((max_size_on_disk != 0) and (local_storage.free_space < size_of_files)):
            DiskSpaceEvent.set()
            if storage.free_space - disk_space_margin < size_of_files:
                logging.warning(
                    "Disk space running low! Starting wait loop ({} s)".format(
                        disk_refresh_time
                    )
                )
                logging.warning(
                    "Disk space: {} GB, Disk space to maintain: {} GiB, ".format(
                        str(float(storage.free_space) / GiB),
                        str(disk_space_margin / GiB),
                    )
                    + f"Restore size: {str(size_of_files / GiB)} GiB"
                )
            elif local_storage.free_space < size_of_files:
                logging.warning(
                    "Reached size limit of files on disk ({} GiB).".format(
                        str(float(max_size_on_disk) / GiB)
                    )
                    + f" Starting wait loop ({disk_refresh_time} s)"
                )
            while (storage.free_space - disk_space_margin < size_of_files) or (
                local_storage.free_space < size_of_files
            ):
                time.sleep(disk_refresh_time)
                # If we can't get disk space, disregard it in the process
                try:
                    storage = StorageState(restore_path)
                except DiskSpaceError:
                    storage = StorageState(restore_path, 0, 0, 0)
                    disk_space_margin = 0
                local_storage = local_storage_for_rstdmf(
                    files_for_rstdmf, restore_path, max_size_on_disk
                )
        DiskSpaceEvent.clear()
        try:
            rstdmf_rename(list_of_path, restore_path)
        except RstdmfError:
            raise


class RstdmfError(Exception):
    pass
