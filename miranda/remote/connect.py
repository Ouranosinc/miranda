"""Remote Connection Operations module."""
from __future__ import annotations

import logging.config
import warnings
from getpass import getpass
from pathlib import Path

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
    """Connection contextualise class."""

    def __init__(
        self,
        username: str | Path = None,
        host: str | Path = None,
        protocol: str = "sftp",
        *args,
        **kwargs,
    ):
        self.user = username or input("Enter username: ")
        self.host = host or input("Enter host URL: ")
        self._args = list(*args)
        self._kwargs = {**kwargs}
        self.__c = None

        if protocol.lower() in ["sftp", "scp"]:
            self.protocol = protocol.lower()
        else:
            raise ValueError('Protocol must be "sftp" or "scp".')

    def update(self, **kwargs):
        """Update connection keyword arguments.

        Warnings
        --------
        Credentials are not encrypted.
        """
        self._kwargs = kwargs

    def __call__(self, **kwargs):
        """Update keyword arguments on call."""
        self.update(**kwargs)
        return self

    def __str__(self):
        return f"Connection to {self.host} as {self.user}"

    def __repr__(self):
        return f"<{self.__class__.__module__}.{self.__class__.__name__} object at {hex(id(self))}>"

    def connect(self, **kwargs):
        """Connect to a remote server with credential prompts."""
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
        except Exception as e:
            raise e

    def __enter__(self, **kwargs):
        return self.connect()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__c.close()
