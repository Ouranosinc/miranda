from __future__ import annotations

import numpy as np
import xarray as xr


def test_test_timeseries(timeseries):
    """Test the timeseries fixture."""
    tt = timeseries(np.array([3, 4, 20, 20, 0, 6, 9, 25, 0, 0]), variable="tas")

    assert tt.name == "tas"
    assert tt.attrs["description"] == "Mean surface temperature."
    assert tt.attrs["standard_name"] == "air_temperature"
    assert tt.attrs["cell_methods"] == "time: mean"
    assert tt.attrs["units"] == "K"
    assert tt.time.size == 10

    resample = tt.resample(time="7D").mean(dim="time")
    np.testing.assert_almost_equal(resample.values, ([8.85714286, 8.33333333]))


def test_multivariable_dataset(multivariable_dataset):
    """Test the multivariable_timeseries fixture."""
    variables = [
        "tasmax",
        "tasmin",
        "tdps",
        "evspsblpot",
        "pr",
        "rdrs",
        "hurs",
        "ua",
        "va",
    ]

    values = np.random.rand(10)
    ds = multivariable_dataset(values, variables=variables)

    assert isinstance(ds, xr.Dataset)
    assert all(var in ds for var in variables)
    assert ds.attrs["title"] == "Multivariable Test Time Series"
    assert ds.attrs["description"] == "This is a test time series with multiple variables."
