#!/usr/bin/env python
"""Tests for `miranda` package."""
from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path

import miranda


class TestMirandaVersion:
    def test_version(self):
        assert isinstance(miranda.__version__, str)


class TestMiranda:
    @classmethod
    def setup_class(cls):
        pass

    def test_something(self):
        pass

    @classmethod
    def teardown_class(cls):
        pass


# class TestDatabase:
#     def test_create_database(self):
#         common = Path(__file__).parent
#         db = miranda.DataBase(common)
#
#         assert len(db) == 3
#         assert str(db.__dict__["_common_path"]).endswith("tests/data/cmip5")
#
#     def test_dict_funcs(self):
#         common = Path(__file__).parent
#         db = miranda.DataBase(common)
#
#         true_keys = set(db.__dict__.keys())
#         assert {"_files", "_is_server", "_source", "_destination"}.issubset(true_keys)
#
#         keys = set(db.keys())
#         assert {
#             "project_name",
#             "recursive",
#             "successful_transfers",
#             "file_suffixes",
#         }.issubset(keys)
#         assert not {"_files", "_is_server", "_source", "_destination"}.issubset(keys)
#
#     def test_url_validator(self):
#         common = Path(__file__).parent
#         db = miranda.DataBase(common)
#
#         url = "https://www.google.ca"
#         short_url = "http://bit.ly/1a2b3c4d5e"
#         not_url = "htttp://not-a-url.biz"
#         assert db._url_validate(url)
#         assert db._url_validate(short_url)
#         assert not db._url_validate(not_url)


def test_package_metadata():
    """Test the package metadata."""
    project = find_spec("miranda").submodule_search_locations[0]

    metadata = Path(project).resolve().joinpath("__init__.py")

    with metadata.open() as f:
        contents = f.read()
        assert """Trevor James Smith""" in contents
        assert '__email__ = "smith.trevorj@ouranos.ca"' in contents
        assert '__version__ = "0.6.0-dev.8"' in contents
