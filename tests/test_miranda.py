from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path

import miranda


class TestMirandaVersion:
    def test_version(self):
        assert isinstance(miranda.__version__, str)


def test_package_metadata():
    """Test the package metadata."""
    project = find_spec("miranda").submodule_search_locations[0]

    metadata = Path(project).resolve().joinpath("__init__.py")

    with metadata.open() as f:
        contents = f.read()
        assert """Trevor James Smith""" in contents
        assert '__email__ = "smith.trevorj@ouranos.ca"' in contents
        assert '__version__ = "0.6.0-dev.20"' in contents
