"""Database Management module."""
from __future__ import annotations

import logging.config
import os
from pathlib import Path
from types import GeneratorType

from .io import find_filepaths
from .scripting import LOGGING_CONFIG
from .units import GiB
from .validators import url_validate

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["DataBase"]


class DataBase:
    """Database management class."""

    def __init__(
        self,
        source,
        *,
        destination: Path | str | None = None,
        common_path: Path | str | None = None,
        file_pattern: str | list[str] = "*.nc",
        project_name: str = None,
        recursive: bool = True,
    ):
        self._source = Path(source)

        if destination is not None:
            self._destination = Path(destination)
        else:
            self._destination = Path().cwd()

        self.project_name = str(project_name)
        if not self.project_name:
            self.project_name = self._destination.stem

        if not file_pattern:
            self.file_suffixes = ["*"]

        elif isinstance(file_pattern, str):
            self.file_suffixes = [file_pattern]
        elif isinstance(file_pattern, (GeneratorType, list)):
            self.file_suffixes = file_pattern

        if not recursive:
            self.recursive = False
        else:
            self.recursive = True

        # if common_path is None:
        #     self._common_path = Path(source)

        self._files = self._scrape(source)
        self._is_server = False

        self.successful_transfers = int(0)

    def __repr__(self):
        """Repl function."""
        return "<{}.{} object at {}>".format(
            self.__class__.__module__, self.__class__.__name__, hex(id(self))
        )

    def __str__(self):
        """String function."""
        prepr = "[%s]" % ", ".join([f'{k}: "{v}"' for k, v in self.__dict__.items()])
        return f"{self.__class__.__name__}({prepr})"

    def __getitem__(self, key):
        """Getter."""
        return self.__dict__[key]

    def __setitem__(self, key, value):
        """Setter."""
        self.__dict__[key] = value

    def __delitem__(self, key):
        """Delete item."""
        del self.__dict__[key]

    def __contains__(self, key):
        """Contains function."""
        return key in self.__dict__

    def __len__(self):
        """Length."""
        return len(self._files)

    def _scrape(self, source) -> list[Path]:
        if source is None:
            raise ValueError("No source provided.")
        if isinstance(source, (GeneratorType, list, tuple, str, Path)):
            files = find_filepaths(source, **self._as_dict())
            common_path = os.path.commonpath(files)
            self._files = files
            self._common_path = common_path
            return files
        raise ValueError("Source must be an iterable of strings or Paths.")

    def _as_dict(self):
        return {
            key: value
            for key, value in self.__dict__.items()
            if not key.startswith("_") and not callable(key)
        }

    def items(self):
        """Show items."""
        return self._as_dict().items()

    def keys(self):
        """Show keys."""
        return self._as_dict().keys()

    def values(self):
        """Show values."""
        return self._as_dict().values()

    def group_by(
        self,
        *,
        common_path: Path | str = None,
        subdirectories: bool = True,
        dates: bool = True,
        size: int = 10 * GiB,
    ):
        """Grouping meta-function.

        Notes
        -----
        Not yet implemented.

        """
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

    def target(self, target: Path | str):
        """Target directory or server address."""
        self._destination = target
        self._is_server = self._url_validate(target=target)

    @staticmethod
    def _url_validate(target):
        return url_validate(target=target)

    def archive(self):
        """Not yet implemented."""
        pass

    def transfer(self):
        """Not yet implemented."""
        pass
