"""Remote Operations module."""

# FIXME: This module should be moved to its own package for licensing reasons.

from __future__ import annotations

import logging.config
import os
import tarfile
import tempfile
import time
import warnings
from pathlib import Path

import miranda.remote
from miranda.scripting import LOGGING_CONFIG

try:
    import fabric  # noqa
    from paramiko import SFTPClient, SSHClient, SSHException  # noqa
    from scp import SCPClient, SCPException  # noqa
except ImportError:
    warnings.warn(
        f"{__name__} functions require additional dependencies. Please install them with `pip install miranda[remote]`."
    )


logging.config.dictConfig(LOGGING_CONFIG)
__all__ = ["create_archive", "create_remote_directory", "transfer_file"]


def create_remote_directory(
    directory: str | os.PathLike[str] | Path,
    transport: SSHClient | fabric.Connection | miranda.remote.Connection,
) -> None:
    """
    Call "mkdir -p" function to create a folder structure over SFTP/SSH and wait for confirmation before continuing.

    Parameters
    ----------
    directory : str or os.PathLike or Path
        The directory to create.
    transport : SSHClient or fabric.Connection or miranda.remote.Connection
        The transport to use.

    Raises
    ------
    ConnectionError
        When the transport is not a valid connection.
    """
    if isinstance(directory, str):
        directory = Path(directory)
    msg = f"Creating remote path: {directory}."
    logging.info(msg)

    ownership = "0775"
    command = f"mkdir -p -m {ownership} '{directory.as_posix()}'"
    if isinstance(transport, (fabric.Connection, miranda.remote.Connection)):
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
    source_files: list[str | os.PathLike[str] | Path],
    destination: str | os.PathLike[str],
    transport: (
        SCPClient | SFTPClient | fabric.Connection | miranda.remote.Connection | None
    ) = None,
    delete: bool = True,
    compression: bool = False,
    recursive: bool = True,
) -> None:
    """
    Create an archive from source files and transfer to another location (remote or local).

    Parameters
    ----------
    source_files : list of str or os.PathLike
        The source files to archive.
    destination : str or os.PathLike
        The destination directory to save the archive.
    transport : SCPClient or SFTPClient or fabric.Connection or miranda.remote.Connection, optional
        The transport to use.
    delete : bool
        Whether to delete the temporary file. Default: True.
    compression : bool
        Whether to compress the archive. Default: False.
    recursive : bool
        Whether to search for files recursively. Default: True.

    Raises
    ------
    ValueError
        If the compression value is invalid.
    """
    if compression:
        write = "w:gz"
    elif not compression:
        write = "w"
    else:
        raise ValueError(f"Compression: {compression}")

    with tempfile.NamedTemporaryFile(delete=delete) as temp:
        archive_file = temp.name
        with tarfile.open(archive_file, write) as tar:
            for name in source_files:
                try:
                    msg = f"Tarring {Path(name).name}."
                    logging.info(msg)
                    tar.add(Path(name).relative_to(Path.cwd()), recursive=recursive)
                except OSError as e:  # noqa: PERF203
                    msg = f'File "{Path(name).name}" failed to be tarred: {e}'
                    logging.warning(msg)
            tar.close()
        transfer_file(archive_file, destination, transport)
    return


def transfer_file(
    source_file: str | os.PathLike[str] | Path,
    destination_file: str | os.PathLike[str] | Path,
    transport: (
        SCPClient | SFTPClient | fabric.Connection | miranda.remote.Connection | None
    ) = None,
) -> bool:
    """
    Transfer file from one location (remote or local) to another.

    Parameters
    ----------
    source_file : str or os.PathLike or Path
        The source file to transfer.
    destination_file : str or os.PathLike or Path
        The destination file to transfer to.
    transport : SCPClient or SFTPClient or fabric.Connection or miranda.remote.Connection, optional
        The transport to use.

    Returns
    -------
    bool
        Whether the transfer was successful.

    Raises
    ------
    SCPException
        If the SCP transfer fails.
    SSHException
        If the SSH connection fails.
    """
    source_file = Path(source_file)
    destination_file = Path(destination_file)

    if transport:
        try:
            msg = f"Beginning transfer of {source_file}."
            logging.info(msg)

            transport.put(str(source_file), str(destination_file))

            msg = f"Transferred { Path(destination_file).name} to {Path(destination_file).parent}."
            logging.info(msg)

        except (OSError, SCPException, SSHException) as e:
            msg = f'File "{destination_file.name}" failed to be transferred: {e}.'
            logging.warning(msg)
            return False

        msg = f"Transferred {Path(destination_file).name} to {Path(destination_file).parent}."
        logging.info(msg)

    else:
        try:
            destination_file.write_bytes(source_file.read_bytes())
        except (OSError, SCPException, SSHException) as e:
            msg = f'File "{source_file.name}" failed to be copied: {e}'
            logging.error(msg)
            return False
    return True
