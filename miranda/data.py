#!/bin/env python3
"""
Copyright 2019 Trevor James Smith

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import os
import re
from pathlib import Path
from types import GeneratorType
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

from .utils import find_filepaths
from .utils import GiB

# from .archive import *
# from .server import Connection

__all__ = ["DataBase"]


class DataBase(object):
    """
    """

    def __init__(
        self,
        source,
        *,
        destination: Optional[Union[Path, str]] = None,
        common_path: Optional[Union[Path, str]] = None,
        file_pattern: str = "*.nc",
        project_name: str = None,
        recursive: bool = True
    ):
        self.destination = Path(destination)
        if not self.destination:
            self.destination = Path().cwd()

        self.project_name = str(project_name)
        if not self.project_name:
            self.project_name = str(self.destination.name)

        self.file_suffixes = str(file_pattern)
        self.recursive = recursive

        self.common_path = Path(source)
        self._files = self._scrape(source)
        if common_path:
            self.common_path = Path(common_path)
        self._is_server = False
        self.source = Path(source)
        self.successful_transfers = int(0)

    def __repr__(self):
        return "<{}.{} object at {}>".format(
            self.__class__.__module__, self.__class__.__name__, hex(id(self))
        )

    def __str__(self):
        prepr = "[%s]" % ", ".join(['{}: "{}"'.format(k, v) for k, v in self.items()])
        return "{}({})".format(self.__class__.__name__, prepr)

    def _scrape(self, source) -> List[Path]:
        if source is None:
            raise ValueError("Source must be a string or Path.")
        elif isinstance(source, (GeneratorType, List, Tuple, str, Path)):
            files = find_filepaths(source, **self._as_dict())
            common_path = os.path.commonpath(f for f in files)
            self._files = files
            self.common_path.update(common_path)
            return files
        else:
            raise ValueError

    def _as_dict(self):
        return {
            key: value
            for key, value in self.__dict__.items()
            if not key.startswith("__") and not callable(key)
        }

    def items(self):
        return self._as_dict().items()

    def keys(self):
        return self._as_dict().keys()

    def values(self):
        return self._as_dict().values()

    def group_by(
        self,
        *,
        common_path: Union[Path, str] = None,
        subdirectories: bool = True,
        dates: bool = True,
        size: int = 10 * GiB
    ):
        # use_grouping = True
        #
        # if subdirectories:
        #     file_groups = group_by_subdirectories(self._files, within=common_path)
        #
        # else:
        #     file_groups = defaultdict(lambda: list())
        #     for f in self._files:
        #         file_groups["."].append(f)
        pass

    def target(self, target: Union[Path, str]):
        self.destination = target
        self._is_server = self._url_validate(target=target)

    @staticmethod
    def _url_validate(target):
        """
        see: https://stackoverflow.com/a/7160778/7322852
        """
        regex = re.compile(
            r"^(?:http|ftp)s?://"  # http:// or https://
            # domain...
            r"(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|"
            r"localhost|"  # localhost...
            r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"  # ...or ip
            r"(?::\d+)?"  # optional port
            r"(?:/?|[/?]\S+)$",
            re.IGNORECASE,
        )
        return re.match(regex, target)

    def archive(self):
        pass

    def transfer(self):
        pass
