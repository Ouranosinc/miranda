from getpass import getpass
from pathlib import Path
from typing import Union

import fabric
from paramiko import SSHClient
from scp import SCPClient

__all__ = ["Connection"]


class Connection(object):
    def __init__(
        self,
        username: Union[str, Path] = None,
        server: Union[str, Path] = None,
        protocol: str = "sftp",
    ):
        self.user = username or input("Enter username: ")
        self.server = server or input("Enter server URL: ")
        self.__c = None

        if protocol.lower() in ["sftp", "scp"]:
            self.protocol = protocol.lower()
        else:
            raise ValueError("Protocol must be 'sftp' or 'scp'.")

    def __str__(self):
        return "Connection to {} as {}".format(self.server, self.user)

    def __repr__(self):
        return "<{}.{} object at {}>".format(
            self.__class__.__module__, self.__class__.__name__, hex(id(self))
        )

    def connect(self):
        try:
            if self.protocol == "sftp":
                c = fabric.Connection(
                    host=self.server,
                    user=self.user,
                    connect_kwargs=dict(password=getpass("Enter password: ")),
                )
                self.__c = c
            else:
                c = SSHClient()
                c.connect(
                    self.server,
                    username=self.user,
                    password=getpass("Enter password: "),
                )
                self.__c = SCPClient(c.get_transport())
            return self.__c
        except Exception as e:
            raise e

    def __enter__(self):
        return self.connect()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__c.close()
