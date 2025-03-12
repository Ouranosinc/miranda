"""Remote Connection Operations module."""

# FIXME: This module should be moved to its own package for licensing reasons.

from __future__ import annotations

import logging.config
import warnings
from getpass import getpass
from pathlib import Path
from typing import Any

from miranda.scripting import LOGGING_CONFIG

try:
    import fabric  # noqa
    from paramiko import SSHClient  # noqa
    from scp import SCPClient  # noqa
except ImportError:
    warnings.warn(
        f"{__name__} functions require additional dependencies."
        f"Please install them with `pip install miranda[remote]`."
    )


logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["Connection"]


class Connection:
    r"""
    Connection contextualise class.

    Parameters
    ----------
    username : str or Path, optional
        The username to use for the connection.
    host : str or Path, optional
        The host URL to connect to.
    protocol : str, optional
        The protocol to use for the connection.
    \*args : list
        Additional arguments.
    \*\*kwargs : dict
        Additional keyword arguments.

    Raises
    ------
    ValueError
        When the protocol is not "sftp" or "scp".

    Warnings
    --------
    Credentials are not encrypted.
    """

    def __init__(
        self,
        username: str | Path | None = None,
        host: str | Path | None = None,
        protocol: str = "sftp",
        *args,
        **kwargs,
    ):
        r"""
        Initialise the connection object.

        Parameters
        ----------
        username : str or Path, optional
            The username to use for the connection.
        host : str or Path, optional
            The host URL to connect to.
        protocol : str, optional
            The protocol to use for the connection.
        \*args : list
            Additional arguments.
        \*\*kwargs : dict
            Additional keyword arguments.

        Raises
        ------
        ValueError
            When the protocol is not "sftp" or "scp".
        """
        self.user = username or input("Enter username: ")
        self.host = host or input("Enter host URL: ")
        self._args = list(*args)
        self._kwargs = {**kwargs}
        self.__c = None

        if protocol.lower() in ["sftp", "scp"]:
            self.protocol = protocol.lower()
        else:
            raise ValueError('Protocol must be "sftp" or "scp".')

    def update(self, **kwargs: dict[str, Any]):
        r"""
        Update connection keyword arguments.

        Parameters
        ----------
        \*\*kwargs : dict
            The keyword arguments to update.

        Warnings
        --------
        Credentials are not encrypted.
        """
        self._kwargs = kwargs

    def __call__(self, **kwargs: dict[str, Any]):
        r"""
        Update keyword arguments on call.

        Parameters
        ----------
        \*\*kwargs : dict
            The keyword arguments to update.

        Returns
        -------
        Connection
            The updated connection object.
        """
        self.update(**kwargs)
        return self

    def __str__(self):
        """
        The string representation of the connection.

        Returns
        -------
        str
            The connection string.
        """
        return f"Connection to {self.host} as {self.user}"

    def __repr__(self):  # noqa: D105
        return f"<{self.__class__.__module__}.{self.__class__.__name__} object at {hex(id(self))}>"

    def connect(self, **kwargs: dict[str, Any]):
        r"""
        Connect to a remote server with credential prompts.

        Parameters
        ----------
        \*\*kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        fabric.Connection or SCPClient
            The connection object.

        Raises
        ------
        Exception
            If the connection fails.
        """
        try:
            keywords = (
                dict(**kwargs)
                or dict(**self._kwargs)
                or dict(password=getpass("Enter password: "))
            )
            if self.protocol == "sftp":
                c = fabric.Connection(
                    host=self.host, user=self.user, connect_kwargs=keywords
                )
                self.__c = c
            else:
                c = SSHClient()
                c.connect(
                    self.host,
                    username=self.user,
                    password=self._kwargs["password"] or getpass("Enter password: "),
                )
                self.__c = SCPClient(c.get_transport())

            return self.__c
        # FIXME: This is too broad.
        except Exception as e:
            raise e

    def __enter__(self, **kwargs):  # noqa: D105
        return self.connect()

    def __exit__(self, exc_type, exc_val, exc_tb):  # noqa: D105
        self.__c.close()
