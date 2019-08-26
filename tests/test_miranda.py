"""
Tests for `miranda` module.
"""
import pytest

import miranda
from miranda import __version__


class TestMirandaVersion(object):
    def test_version(self):
        assert __version__


class TestMiranda(object):
    @classmethod
    def setup_class(cls):
        pass

    def test_something(self):
        pass

    @classmethod
    def teardown_class(cls):
        pass
