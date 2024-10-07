"""Archive Module."""

from __future__ import annotations

import logging.config
import os
from collections import defaultdict
from pathlib import Path

from miranda.archive import (
    group_by_deciphered_date,
    group_by_size,
    group_by_subdirectories,
)
from miranda.io import find_filepaths
from miranda.scripting import LOGGING_CONFIG
from miranda.storage import report_file_size
from miranda.utils import single_item_list, working_directory

from .connect import Connection
from .ops import create_archive, create_remote_directory, transfer_file

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["archive_database"]


def archive_database(
    source: Path | str | os.PathLike[str] | list[str | os.PathLike[str] | Path],
    common_path: Path | str | os.PathLike[str],
    destination: Path | str | os.PathLike[str],
    file_suffixes: str = ".nc",
    server: str | None = None,
    username: str | None = None,
    project_name: str | None = None,
    overwrite: bool = False,
    compression: bool = False,
    recursive: bool = False,
    use_grouping: bool = True,
    use_subdirectories: bool = True,
    dry_run: bool = False,
) -> None:
    """
    Archive database files to a remote server.

    Given a source, destination, and dependent on file size limit, create tarfile archives and transfer
    files to another server for backup purposes.

    Parameters
    ----------
    source : Path or str or os.PathLike or list
        The source directory containing the files to archive.
    common_path : Path or str or os.PathLike
        The common path to use for grouping files.
    destination : Path or str or os.PathLike
        The destination directory to save the files.
    file_suffixes : str
        The file suffix to use for filtering files.
    server : str, optional
        The server to connect to.
    username : str, optional
        The username to use for the connection.
    project_name : str, optional
        The project name to use for the files.
    overwrite : bool, optional
        Whether to overwrite existing files.
    compression : bool, optional
        Whether to compress the files.
    recursive : bool, optional
        Whether to search for files recursively.
    use_grouping : bool, optional
        Whether to group files by date.
    use_subdirectories : bool, optional
        Whether to use subdirectories for grouping.
    dry_run : bool, optional
        Whether to run in dry-run mode.

    Raises
    ------
    RuntimeError
        If the transfer fails.
    """
    project = "{project_name}_{group_name}_{common_date}{part}.{suffix}"

    if not project_name:
        project_name = destination.name

    if compression:
        suffix = "tar.gz"
    elif not compression:
        suffix = "tar"
    else:
        raise ValueError(f"Compression: {compression}")

    file_list, source_path = find_filepaths(
        source=source, recursive=recursive, file_suffixes=file_suffixes
    )

    if use_subdirectories:
        file_groups = group_by_subdirectories(file_list, within=common_path)

    else:
        file_groups = defaultdict(list)
        for f in file_list:
            file_groups["."].append(f)

    connection = Connection(username=username, host=server, protocol="sftp")

    if dry_run:
        logging.info(
            "Running archival functions in `dry_run` mode. No files will be transferred."
        )

    try:
        successful_transfers = list()
        with connection as ctx:
            for group_name, members in file_groups.items():
                remote_path = Path(destination, group_name)

                if not dry_run:
                    if not remote_path.exists():
                        create_remote_directory(remote_path, transport=ctx)

                if use_grouping:
                    dated_groups = group_by_deciphered_date(members)
                else:
                    dated_groups = dict(group_name=members)

                for common_date, files in dated_groups.items():
                    if not use_grouping or single_item_list(files):
                        for archive_file in files:
                            transfer = Path(remote_path, archive_file.name)

                            if transfer.is_file():
                                if not overwrite:
                                    msg = f"{transfer} exists. Skipping file."
                                    logging.info(msg)
                                    continue
                                msg = f"{transfer} exists. Overwriting."
                                logging.info(msg)
                            if not dry_run:
                                if transfer_file(archive_file, transfer, transport=ctx):
                                    successful_transfers.append(archive_file)
                            else:
                                successful_transfers.append(archive_file)

                    elif use_grouping or not single_item_list(files):
                        sized_groups = group_by_size(files)

                        for i, sized_group in enumerate(sized_groups):
                            if len(sized_groups) > 1:
                                part = f"_{str(i + 1).zfill(3)}"
                            else:
                                part = ""

                            archive_file = project.format(
                                project_name, group_name, common_date, part, suffix
                            )
                            transfer = Path(remote_path, archive_file)

                            if transfer.is_file():
                                if not overwrite:
                                    msg = f'File "{transfer}" exists. Skipping file.'
                                    logging.info(msg)
                                    continue
                                msg = f'File "{transfer}" exists. Overwriting.'
                                logging.info(msg)

                            with working_directory(source_path):
                                if not dry_run:
                                    if create_archive(
                                        sized_group,
                                        transfer,
                                        transport=ctx,
                                        compression=compression,
                                        recursive=recursive,
                                    ):
                                        successful_transfers.extend(sized_group)
                                else:
                                    successful_transfers.extend(sized_group)
                    else:
                        raise FileNotFoundError("No files found in grouping.")

        msg = (
            f"Transferred {len(successful_transfers)} "
            f"of {len([f for f in file_list])} files "
            f"totalling {report_file_size(successful_transfers)}."
        )
        logging.info(msg)

    except Exception as e:
        msg = f"{e}: Failed to transfer files."
        logging.error(msg)
        raise RuntimeError(msg) from e
