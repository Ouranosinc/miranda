from __future__ import annotations
import pathlib
from importlib.util import find_spec

import miranda


class TestMirandaVersion:
    def test_version(self):
        assert isinstance(miranda.__version__, str)


def test_package_metadata():
    """Test the package metadata."""
    project = find_spec("miranda")

    assert project is not None
    assert project.submodule_search_locations is not None
    location = project.submodule_search_locations[0]

    metadata = pathlib.Path(location).resolve().joinpath("__init__.py")

    with metadata.open() as f:
        contents = f.read()
        assert """Trevor James Smith""" in contents
        assert '__email__ = "smith.trevorj@ouranos.ca"' in contents
        assert '__version__ = "0.6.0"' in contents
