from __future__ import annotations

import re

import pytest  # noqa

try:
    import fabric  # noqa
except ImportError:
    fabric = None

from miranda.remote.connect import Connection


@pytest.mark.skipif(not fabric, reason="Needs fabric")
class TestConnection:
    def test_connection_dict(self):
        self.c = Connection(username="qwerty", host="localhost")
        assert ("user", "qwerty") in self.c.__dict__.items()
        assert "localhost" in self.c.__dict__.values()
        assert "protocol" in self.c.__dict__.keys()

    def test_connection_location(self):
        self.c = Connection(username="qwerty", host="localhost")
        match = re.search(r"0x[\d|a-f]+>", repr(self.c))
        assert match is not None

    def test_connection_update(self):
        self.c = Connection(username="qwerty", host="localhost")
        self.c.update(password="12345")
        with self.c as context:
            assert context.is_connected is False

    def test_connection_context_pass(self):
        self.c = Connection(username="qwerty", host="localhost")
        with self.c(password="12345") as context:
            assert context.is_connected is False

    def test_connection(self):
        self.c = Connection(username="qwerty", host="localhost")
        c = self.c.connect(password="12345")
        assert c.is_connected is False
