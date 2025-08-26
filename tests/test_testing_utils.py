from __future__ import annotations

import pytest


@pytest.mark.parametrize(
    "output_format,expected_start,expected_hyperlink",
    [
        (
            "md",
            ".. :changelog:\n\n# Changelog\n\n",
            "[@Zeitsperre](https://github.com/Zeitsperre)",
        ),
        (
            "rst",
            ".. :changelog:\n\n=========\nChangelog\n=========\n\n",
            "`@Zeitsperre <https://github.com/Zeitsperre>`_",
        ),
        ("latex", "NotImplemented", "NotImplemented"),
    ],
    ids=["Markdown", "reStructuredText", "NotImplemented"],
)
def test_publish_release_notes(output_format, expected_start, expected_hyperlink):
    from miranda.testing.utils import publish_release_notes

    if expected_start == "NotImplemented":
        with pytest.raises(NotImplementedError):
            publish_release_notes(output_format)
    else:
        release_notes = publish_release_notes(output_format)
        assert isinstance(release_notes, str)
        assert expected_start in release_notes
        assert expected_hyperlink in release_notes


def test_show_versions():
    from miranda.testing.utils import show_versions

    versions_found = show_versions()
    min_dependencies = [
        "cftime",
        "cf-xarray",
        "dask",
        "distributed",
        "h5netcdf",
        "netCDF4",
        "numcodecs",
        "numpy",
        "pandas",
        "pyessv",
        "pyyaml",
        "schema",
        "xarray",
        "xclim",
        "zarr",
    ]

    assert isinstance(versions_found, str)
    missing_deps = [library for library in min_dependencies if f"{library} : None" in versions_found]
    if missing_deps:
        raise ImportError(f"Missing dependencies: {', '.join(missing_deps)}. Try installing `miranda` or adjusting deps.")
