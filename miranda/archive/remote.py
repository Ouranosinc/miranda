import logging
import tarfile
import tempfile
import time
from pathlib import Path
from typing import List
from typing import Union

import fabric
from paramiko import SFTPClient
from paramiko import SSHClient
from paramiko import SSHException
from scp import SCPClient
from scp import SCPException

from miranda.server import Connection

# from functools import singledispatch

__all__ = ["create_archive", "create_remote_directory", "transfer_file"]


def create_remote_directory(
    directory, transport: Union[SSHClient, fabric.Connection, Connection]
):
    """
    This calls a function to create a folder structure over SFTP/SSH and waits
     for confirmation before continuing
    """
    logging.info("Creating remote path: {}".format(directory))

    ownership = "0775"
    command = "mkdir -p -m {} '{}'".format(ownership, directory)
    if isinstance(transport, (fabric.Connection, Connection)):
        with transport:
            transport.run(command)
    elif isinstance(transport, (SSHClient, SCPClient)):
        transport.exec_command(command, timeout=1)
        for i in range(5):
            if not directory.exists():
                time.sleep(1)
                continue
            break
    else:
        raise ConnectionError
    return


def create_archive(
    source_files: List[Union[Path, str]],
    destination: Union[Path, str],
    transport: Union[SCPClient, SFTPClient, fabric.Connection, Connection] = None,
    delete: bool = True,
    compression: bool = False,
    recursive: bool = True,
) -> None:
    if compression:
        write = "w:gz"
    elif not compression:
        write = "w"
    else:
        raise ValueError("Compression: {}".format(compression))

    with tempfile.NamedTemporaryFile(delete=delete) as temp:
        archive_file = temp.name
        with tarfile.open(archive_file, write) as tar:
            for name in source_files:
                try:
                    logging.info("Tarring {}".format(name.name))
                    tar.add(name.relative_to(Path.cwd()), recursive=recursive)
                except Exception as e:
                    msg = 'File "{}" failed to be tarred: {}'.format(name.name, e)
                    logging.warning(msg)
            tar.close()
        transfer_file(archive_file, destination, transport)
    return


def transfer_file(
    source_file: Union[Path, str],
    destination_file: Union[Path, str],
    transport: Union[SCPClient, SFTPClient, fabric.Connection, Connection] = None,
) -> bool:

    source_file = Path(source_file)
    destination_file = Path(destination_file)

    if transport:
        try:
            logging.info("Beginning transfer of {}".format(source_file))
            transport.put(str(source_file), str(destination_file))
            logging.info(
                "Transferred {} to {}".format(
                    Path(destination_file).name, Path(destination_file).parent
                )
            )

        except SCPException or SSHException or IOError or OSError as e:
            msg = 'File "{}" failed to be transferred: {}.'.format(
                destination_file.name, e
            )
            logging.warning(msg)
            return False

        logging.info(
            "Transferred {} to {}".format(
                Path(destination_file).name, Path(destination_file).parent
            )
        )

    else:
        try:
            destination_file.write_bytes(source_file.read_bytes())
        except Exception as e:
            msg = 'File "{}" failed to be copied: {}'.format(source_file.name, e)
            logging.error(msg)
            return False
    return True
