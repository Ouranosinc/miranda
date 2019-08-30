#!/bin/env python3
import logging
from collections import defaultdict
from datetime import datetime as dt
from getpass import getpass
from pathlib import Path
from types import GeneratorType
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import fabric
from paramiko import AuthenticationException

from miranda.archive._utils import _find_files
from miranda.archive._utils import _make_remote_directory
from miranda.archive._utils import _transfer_archive
from miranda.archive._utils import _transfer_single
from miranda.archive.generic_archiver import file_size
from miranda.archive.generic_archiver import group_by_deciphered_date
from miranda.archive.generic_archiver import group_by_size
from miranda.archive.generic_archiver import group_by_subdirectories
from miranda.archive.generic_archiver import single_item_list
from miranda.archive.generic_archiver import working_directory


class DataBase(object):
    """

    """

    def __init__(
        self,
        source,
        *,
        destination: Optional[Union[str or Path]] = None,
        common_path: Optional[Union[str or Path]] = None,
        file_pattern: str = "*.nc",
        project_name: str = None,
        recursive: bool = True,
    ):
        self.destination = destination
        if not self.destination:
            self.destination = Path().cwd()

        self.project_name = project_name
        if not self.project_name:
            self.project_name = self.destination.name

        self.file_suffixes = file_pattern
        self.recursive = recursive
        self._files, self.common_path = self._scrape(source)
        if common_path:
            self.common_path = common_path

        self.source = Path(source)

    def __repr__(self):
        return "<{}.{} object at {}>".format(
            self.__class__.__module__, self.__class__.__name__, hex(id(self))
        )

    def __str__(self):
        prepr = "[%s]" % ", ".join(['{}: "{}"'.format(k, v) for k, v in self.items()])
        return "{}({})".format(self.__class__.__name__, prepr)

    def _scrape(self, source):
        if source is None:
            raise ValueError("Source must be a string or Path.")
        elif isinstance(source, (GeneratorType, List, Tuple, str, Path)):
            files, source = _find_files(source, **self._as_dict())
            self._files = files
            return files, source
        else:
            raise ValueError

    def _as_dict(self):
        return {
            key: value
            for key, value in self.__dict__.items()
            if not key.startswith("__") and not callable(key)
        }

    def items(self):
        return self._as_dict().items()

    def keys(self):
        return self._as_dict().keys()

    def values(self):
        return self._as_dict().values()

    def group_by(
        self,
        *,
        common_path: str or Path = None,
        subdirectories: bool = True,
        dates: bool = True,
    ):
        # use_grouping = True
        #
        # if subdirectories:
        #     file_groups = group_by_subdirectories(self._files, within=common_path)
        #
        # else:
        #     file_groups = defaultdict(lambda: list())
        #     for f in self._files:
        #         file_groups["."].append(f)
        pass

    def target(self, target: str or Path):
        self.destination = target

    def archive(self):
        pass


def archive(
    source: str or Path or List,
    common_path: str or Path,
    destination: str or Path,
    file_suffixes: str = ".nc",
    server: str = None,
    username: str = None,
    password: str = None,
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

    file_list, source_path = _find_files(
        source=source, recursive=recursive, file_suffixes=file_suffixes
    )

    if use_subdirectories:
        file_groups = group_by_subdirectories(file_list, within=common_path)

    else:
        file_groups = defaultdict(lambda: list())
        for f in file_list:
            file_groups["."].append(f)

    try:
        user = username or input("Enter username: ")
        pw = password or getpass("Enter password: ")
        connection = fabric.Connection(
            host=server, user=user, connect_kwargs=dict(password=pw)
        )
    except AuthenticationException as e:
        logging.error("{}: Unable to connect to remote host {}.".format(e, server))
        raise

    try:
        successful_transfers = list()
        with connection as ctx:
            for group_name, members in file_groups.items():
                remote_path = Path(destination, group_name)

                if not remote_path.exists():
                    _make_remote_directory(remote_path, transport=ctx)

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
                                        "{}: {} exists. Skipping file.".format(
                                            dt.now().strftime("%Y-%m-%d %X"), transfer
                                        )
                                    )
                                    continue
                                logging.info(
                                    "{}: {} exists. Overwriting.".format(
                                        dt.now().strftime("%Y-%m-%d %X"), transfer
                                    )
                                )

                            if _transfer_single(archive_file, transfer, transport=ctx):
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
                                        "{}: {} exists. Skipping file.".format(
                                            dt.now().strftime("%Y-%m-%d %X"), transfer
                                        )
                                    )
                                    continue
                                logging.info(
                                    "{}: {} exists. Overwriting.".format(
                                        dt.now().strftime("%Y-%m-%d %X"), transfer
                                    )
                                )

                            with working_directory(source_path):
                                if _transfer_archive(
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
            "Transferred {} of {} files totalling {} at {}.".format(
                len(successful_transfers),
                len([f for f in file_list]),
                file_size(successful_transfers),
                dt.now().strftime("%Y-%m-%d %X"),
            )
        )

    except Exception as e:
        msg = "{}: {} Failed to transfer files.".format(
            dt.now().strftime("%Y-%m-%d %X"), e
        )
        logging.error(msg)
        raise RuntimeError(msg) from e

    return
