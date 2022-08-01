import logging
import logging.config
import tarfile
import tempfile
import time
from pathlib import Path
from typing import List, Union

from miranda.scripting import LOGGING_CONFIG

from .connect import Connection

try:
    import fabric  # noqa
    from paramiko import SFTPClient, SSHClient, SSHException  # noqa
    from scp import SCPClient, SCPException  # noqa
except ImportError:
    raise ImportError(
        f"{__name__} functions require additional dependencies. Please install them with `pip install miranda[full]`."
    )


logging.config.dictConfig(LOGGING_CONFIG)
__all__ = ["create_archive", "create_remote_directory", "transfer_file"]


def create_remote_directory(
    directory: Union[str, Path],
    transport: Union[SSHClient, fabric.Connection, Connection],
) -> None:
    """
    This calls a "mkdir -p" function to create a folder structure over SFTP/SSH and waits
    for confirmation before continuing

    Parameters
    ----------
    directory : Union[str, Path]
    transport : Union[SSHClient, fabric.Connection, Connection]

    Returns
    -------
    None

    """
    if isinstance(directory, str):
        directory = Path(directory)
    logging.info(f"Creating remote path: {directory}")

    ownership = "0775"
    command = f"mkdir -p -m {ownership} '{directory.as_posix()}'"
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
    """

    Parameters
    ----------
    source_files: List[Union[Path, str]]
    destination: Union[Path, str]
    transport: Union[SCPClient, SFTPClient, fabric.Connection, Connection]
    delete: bool
    compression: bool = False
    recursive: bool = True

    Returns
    -------
    None

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
                    logging.info(f"Tarring {name.name}")
                    tar.add(name.relative_to(Path.cwd()), recursive=recursive)
                except Exception as e:
                    msg = f'File "{name.name}" failed to be tarred: {e}'
                    logging.warning(msg)
            tar.close()
        transfer_file(archive_file, destination, transport)
    return


def transfer_file(
    source_file: Union[Path, str],
    destination_file: Union[Path, str],
    transport: Union[SCPClient, SFTPClient, fabric.Connection, Connection] = None,
) -> bool:
    """

    Parameters
    ----------
    source_file: Union[Path, str]
    destination_file: Union[Path, str]
    transport: Union[SCPClient, SFTPClient, fabric.Connection, Connection]

    Returns
    -------
    bool

    """

    source_file = Path(source_file)
    destination_file = Path(destination_file)

    if transport:
        try:
            logging.info(f"Beginning transfer of {source_file}")
            transport.put(str(source_file), str(destination_file))
            logging.info(
                f"Transferred { Path(destination_file).name} to {Path(destination_file).parent}"
            )

        except (SCPException, SSHException, OSError) as e:
            msg = f'File "{destination_file.name}" failed to be transferred: {e}.'
            logging.warning(msg)
            return False

        logging.info(
            f"Transferred {Path(destination_file).name} to {Path(destination_file).parent}"
        )

    else:
        try:
            destination_file.write_bytes(source_file.read_bytes())
        except (SCPException, SSHException, OSError) as e:
            msg = f'File "{source_file.name}" failed to be copied: {e}'
            logging.error(msg)
            return False
    return True
