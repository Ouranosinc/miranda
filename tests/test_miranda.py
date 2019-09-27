"""
Tests for `miranda` module.
"""
from pathlib import Path

import miranda
from .common import test_data


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


class TestMakeClasses:
    def test_database(self):
        common = Path(__file__).parent
        db = miranda.DataBase(common)
        # db_keys = db.keys()

        # TODO: Investigate why this doesn't work
        # assert {"_files", "_is_server", "_source", "_destination"}.issubset(db_keys)
        assert len(db.__dict__["_files"]) == 3
        assert len(db.__dict__["_common_path"]) == 1
