from __future__ import annotations
import logging
import os
from collections.abc import Callable
from os.path import commonpath
from pathlib import Path

import pooch
import pytest
from xclim.testing.helpers import test_timeseries

from miranda.testing.utils import (
    TESTDATA_BRANCH,
    TESTDATA_CACHE_DIR,
    TESTDATA_REPO_URL,
    default_testdata_cache,
    gather_testing_data,
    testing_setup_warnings,
)
from miranda.testing.utils import cassini as _cassini
from miranda.testing.utils import open_dataset as _open_dataset


logger = logging.getLogger("miranda")


@pytest.fixture
def timeseries() -> Callable:
    """Fixture to provide a test time series."""
    return test_timeseries


@pytest.fixture
def multivariable_dataset(timeseries) -> Callable:
    """Return a test time series for the specified variable."""

    def _multivariable_dataset(values, variables: list[str], **kwargs):
        """Return a test time series for multiple variables."""
        kwargs.setdefault("as_dataset", True)
        tt = timeseries(values, variable=variables[0], **kwargs)

        del kwargs["as_dataset"]
        for var in variables[1:]:
            tt[var] = timeseries(values, variable=var, **kwargs)
        tt.attrs["title"] = "Multivariable Test Time Series"
        tt.attrs["description"] = "This is a test time series with multiple variables."
        return tt

    return _multivariable_dataset


@pytest.fixture(scope="session")
def threadsafe_data_dir(tmp_path_factory):
    return Path(tmp_path_factory.getbasetemp().joinpath("data"))


@pytest.fixture(scope="session")
def cassini(threadsafe_data_dir, worker_id):
    return _cassini(
        repo=TESTDATA_REPO_URL,
        branch=TESTDATA_BRANCH,
        cache_dir=(TESTDATA_CACHE_DIR if worker_id == "master" else threadsafe_data_dir),
    )


@pytest.fixture(scope="session")
def open_dataset(threadsafe_data_dir, worker_id):
    def _open_session_scoped_file(file: str | os.PathLike, **xr_kwargs):
        cassini_kwargs = {
            "branch": TESTDATA_BRANCH,
            "repo": TESTDATA_REPO_URL,
            "cache_dir": (TESTDATA_CACHE_DIR if worker_id == "master" else threadsafe_data_dir),
        }
        xr_kwargs.setdefault("cache", True)
        xr_kwargs.setdefault("engine", "h5netcdf")
        return _open_dataset(
            file,
            cassini_kwargs=cassini_kwargs,
            **xr_kwargs,
        )

    return _open_session_scoped_file


@pytest.fixture(autouse=True, scope="session")
def gather_session_data(request, cassini, worker_id):
    """
    Gather testing data on pytest run.

    When running pytest with multiple workers, one worker will copy data remotely to the default cache dir while
    other workers wait using a lockfile. Once the lock is released, all workers will then copy data to their local
    threadsafe_data_dir. As this fixture is scoped to the session, it will only run once per pytest run.
    """
    testing_setup_warnings()
    gather_testing_data(worker_cache_dir=cassini.path, worker_id=worker_id)

    def remove_data_written_flag():
        """Clean up the cache folder once we are finished."""
        flag = default_testdata_cache.joinpath(".data_written")
        if flag.exists():
            try:
                flag.unlink()
            except FileNotFoundError:
                logger.info("Teardown race condition occurred: .data_written flag already removed. Lucky!")
                pass

    request.addfinalizer(remove_data_written_flag)


@pytest.fixture(scope="session")
def era5_precip(cassini):
    """ERA5 precipitation data."""
    era5_precip = "ECMWF/era5_tp_ptype_2024.zip"
    test_data_path = cassini.fetch(era5_precip, downloader=pooch.Unzip())
    common = Path(commonpath(test_data_path))

    data = {
        "tp": common / "data_stream-oper_stepType-accum.nc",
        "ptype": common / "data_stream-oper_stepType-instant.nc",
    }

    return data
