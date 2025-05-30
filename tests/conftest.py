from __future__ import annotations

from pathlib import Path
from typing import Callable

import pytest
from xclim.testing.helpers import test_timeseries

data_path = Path.cwd()

test_data = dict()
test_data["current"] = data_path / "cmip5" / "tmax.cur.nc"
test_data["future"] = data_path / "cmip5" / "tmax.fut.nc"
test_data["observed"] = data_path / "cmip5" / "tmax.obs.nc"


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
