import logging
import tarfile
import tempfile
import time
from datetime import datetime as dt
from pathlib import Path
from typing import List
from typing import Union

import fabric
from paramiko import SFTPClient
from paramiko import SSHClient
from paramiko import SSHException
from scp import SCPClient
from scp import SCPException


__all__ = ["_transfer_single", "_transfer_archive", "_make_remote_directory"]


def _make_remote_directory(directory, transport: Union[SSHClient, fabric.Connection]):
    """
    This calls a function to create a folder structure over SFTP/SSH and waits
     for confirmation before continuing
    """
    logging.info(
        "{}: Creating remote path: {}".format(
            dt.now().strftime("%Y-%m-%d %X"), directory
        )
    )

    ownership = "0775"
    command = "mkdir -p -m {} '{}'".format(ownership, directory)
    if isinstance(transport, fabric.Connection):
        transport.run(command)
    elif isinstance(transport, (SSHClient, SCPClient)):
        transport.exec_command(command, timeout=1)
        for i in range(5):
            if not directory.exists():
                time.sleep(1)
                continue
            break
    return


def _transfer_archive(
    source_files: List[Union[Path, str]],
    destination: Union[Path, str],
    transport: Union[SCPClient, SFTPClient, fabric.Connection],
    compression: bool = False,
    recursive: bool = True,
) -> bool:
    if compression:
        write = "w:gz"
    elif not compression:
        write = "w"
    else:
        raise ValueError("Compression: {}".format(compression))

    with tempfile.NamedTemporaryFile() as temp:
        archive_file = temp.name
        with tarfile.open(archive_file, write) as tar:
            for name in source_files:
                try:
                    logging.info(
                        "{}: Tarring {}".format(
                            dt.now().strftime("%Y-%m-%d %X"), name.name
                        )
                    )
                    tar.add(name.relative_to(Path.cwd()), recursive=recursive)
                except Exception as e:
                    msg = '{}: File "{}" failed to be tarred: {}'.format(
                        dt.now().strftime("%Y-%m-%d %X"), name.name, e
                    )
                    logging.warning(msg)

            tar.close()

        try:
            logging.info(
                "{}: Beginning scp transfer of {} to {}".format(
                    dt.now().strftime("%Y-%m-%d %X"),
                    Path(destination).name,
                    Path(destination).parent,
                )
            )
            transport.put(str(archive_file), str(destination))

        except SCPException or SSHException or IOError or OSError as e:
            msg = '{}: File "{}" failed to be added: {}.'.format(
                dt.now().strftime("%Y-%m-%d %X"), destination.name, e
            )
            logging.warning(msg)
            return False

        logging.info(
            "{}: Transferred {} to {}".format(
                dt.now().strftime("%Y-%m-%d %X"),
                Path(destination).name,
                Path(destination).parent,
            )
        )
    return True


def _transfer_single(
    source_file: Union[Path, str],
    destination: Union[Path, str],
    transport: Union[SCPClient, SFTPClient, fabric.Connection],
) -> bool:
    try:
        logging.info(
            "{}: Passing {}".format(dt.now().strftime("%Y-%m-%d %X"), source_file)
        )
        transport.put(str(source_file), str(destination))
        logging.info(
            "{}: Transferred {} to {}".format(
                dt.now().strftime("%Y-%m-%d %X"),
                Path(destination).name,
                Path(destination).parent,
            )
        )
    except SCPException or SSHException or IOError or OSError as e:
        msg = '{}: File "{}" failed to be added: {}.'.format(
            dt.now().strftime("%Y-%m-%d %X"), destination.name, e
        )
        logging.warning(msg)
        return False
    return True
