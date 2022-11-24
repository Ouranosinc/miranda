import logging.config
import warnings
from getpass import getpass
from pathlib import Path
from typing import Union

from miranda.scripting import LOGGING_CONFIG

try:
    import fabric  # noqa
    from paramiko import SSHClient  # noqa
    from scp import SCPClient  # noqa
except ImportError:
    warnings.warn(
        f"{__name__} functions require additional dependencies."
        f"Please install them with `pip install miranda[full]`."
    )


logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["Connection"]


class Connection:
    def __init__(
        self,
        username: Union[str, Path] = None,
        host: Union[str, Path] = None,
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
        self._kwargs = kwargs

    def __call__(self, **kwargs):
        self.update(**kwargs)
        return self

    def __str__(self):
        return f"Connection to {self.host} as {self.user}"

    def __repr__(self):
        return f"<{self.__class__.__module__}.{self.__class__.__name__} object at {hex(id(self))}>"

    def connect(self, **kwargs):
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
