#!/bin/env python3
import logging
from collections import defaultdict
from datetime import datetime as dt
from pathlib import Path
from typing import List
from typing import Union

from .archive.groupings import group_by_deciphered_date
from .archive.groupings import group_by_size
from .archive.groupings import group_by_subdirectories
from .connect import Connection
from .ops import create_archive
from .ops import create_remote_directory
from .ops import transfer_file
from .storage import report_file_size
from .utils import find_filepaths
from .utils import single_item_list
from .utils import working_directory

__all__ = ["archive_database"]


def archive_database(
    source: Union[Path, str, List],
    common_path: Union[Path, str],
    destination: Union[Path, str],
    file_suffixes: str = ".nc",
    server: str = None,
    username: str = None,
    project_name: str = None,
    overwrite: bool = False,
    compression: bool = False,
    recursive: bool = False,
    use_grouping: bool = True,
    use_subdirectories: bool = True,
) -> None:
    """
    Given a source, destination, and dependent on file size limit, create tarfile archives and transfer
     files to another server for backup purposes
    """
    project = "{}_{}{}.{}"

    if not project_name:
        project_name = destination.name

    if compression:
        suffix = "tar.gz"
    elif not compression:
        suffix = "tar"
    else:
        raise ValueError("Compression: {}".format(compression))

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

    try:
        successful_transfers = list()
        with connection as ctx:
            for group_name, members in file_groups.items():
                remote_path = Path(destination, group_name)

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
                                    logging.info(
                                        "{} exists. Skipping file.".format(transfer)
                                    )
                                    continue
                                logging.info("{} exists. Overwriting.".format(transfer))

                            if transfer_file(archive_file, transfer, transport=ctx):
                                successful_transfers.append(archive_file)

                    elif use_grouping or not single_item_list(files):
                        sized_groups = group_by_size(files)

                        for i, sized_group in enumerate(sized_groups):
                            if len(sized_groups) > 1:
                                part = "_{}".format(str(i + 1).zfill(3))
                            else:
                                part = ""

                            archive_file = project.format(
                                project_name, group_name, common_date, part, suffix
                            )
                            transfer = Path(remote_path, archive_file)

                            if transfer.is_file():
                                if not overwrite:
                                    logging.info(
                                        'File "{}" exists. Skipping file.'.format(
                                            transfer
                                        )
                                    )
                                    continue
                                logging.info(
                                    'File "{}" exists. Overwriting.'.format(transfer)
                                )

                            with working_directory(source_path):
                                if create_archive(
                                    sized_group,
                                    transfer,
                                    transport=ctx,
                                    compression=compression,
                                    recursive=recursive,
                                ):
                                    successful_transfers.extend(sized_group)
                    else:
                        raise FileNotFoundError("No files found in grouping.")

        logging.info(
            "Transferred {} of {} files totalling {}.".format(
                len(successful_transfers),
                len([f for f in file_list]),
                report_file_size(successful_transfers),
            )
        )

    except Exception as e:
        msg = "{}: Failed to transfer files.".format(e)
        logging.error(msg)
        raise RuntimeError(msg) from e


if __name__ == "__main__":
    logging.basicConfig(
        filename="{}_{}.log".format(
            dt.strftime(dt.now(), "%Y%m%d"), Path(__name__).stem
        ),
        level=logging.INFO,
    )
