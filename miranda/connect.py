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
        return "Connection to {} as {}".format(self.host, self.user)

    def __repr__(self):
        return "<{}.{} object at {}>".format(
            self.__class__.__module__, self.__class__.__name__, hex(id(self))
        )

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


if __name__ == "__main__":
    a = Connection(username="tttt", host="1234")
    with a.connect(password="stuff") as f:
        print(f.is_connected)
