from getpass import getpass
from pathlib import Path
from typing import Union

import fabric
from paramiko import SSHClient
from scp import SCPClient


class Connection(object):
    def __init__(
        self,
        username: Union[str, Path] = None,
        server: Union[str, Path] = None,
        protocol: str = "sftp",
    ):
        self.user = username or input("Enter username: ")
        self.server = server or input("Enter server URL: ")

        if protocol.lower() in ["sftp", "scp"]:
            self.protocol = protocol.lower()
        else:
            raise ValueError("Protocol must be 'sftp' or 'scp'.")

    def __repr__(self):
        return "<{}.{} object at {}>".format(
            self.__class__.__module__, self.__class__.__name__, hex(id(self))
        )

    def connect(self):
        if self.protocol == "sftp":
            c = fabric.Connection(
                host=self.server,
                user=self.user,
                connect_kwargs=dict(password=getpass("Enter password: ")),
            )
            return c
        else:
            c = SSHClient()
            c.connect(
                self.server, username=self.user, password=getpass("Enter password: ")
            )
            return SCPClient(c.get_transport())
