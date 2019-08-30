"""
Tests for `miranda` module.
"""
import miranda


class TestMirandaVersion(object):
    def test_version(self):
        assert isinstance(miranda.__version__, str)


class TestMiranda(object):
    @classmethod
    def setup_class(cls):
        pass

    def test_something(self):
        pass

    @classmethod
    def teardown_class(cls):
        pass
