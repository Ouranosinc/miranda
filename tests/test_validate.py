from __future__ import annotations

import pytest
import xarray as xr
from schema import SchemaError

from miranda.validate import cf_dimensions_schema, url_validate


def test_valid_dimensions_schema():
    valid_input = {
        "dimensions": {
            "time": {
                "_cf_dimension_name": "time",
                "_ensure_correct_time": "1YS",
                "_strict_time": True,
                "axis": "T",
                "bounds": "time_bnds",
                "standard_name": "time",
                "long_name": "Time",
                "units": "days since 1850-01-01",
            },
            "lat": {
                "_cf_dimension_name": "lat",
                "_precision": 2,
                "axis": "Y",
                "standard_name": "latitude",
                "long_name": "Latitude",
                "units": "degrees_north",
            },
            "lon": {
                "_cf_dimension_name": "lon",
                "_precision": 2,
                "axis": "X",
                "standard_name": "longitude",
                "long_name": "Longitude",
                "units": "degrees_east",
            },
        }
    }
    assert cf_dimensions_schema.validate(valid_input)


def test_invalid_time_schema():
    invalid_time = {
        "dimensions": {
            "time": {
                "_cf_dimension_name": "time",
                "_ensure_correct_time": "1D",
                "_strict_time": False,
                "axis": "T",
                "standard_name": "time",
                "long_name": "Tijd",
                "units": "days since 2000-01-01",
            }
        }
    }
    try:
        cf_dimensions_schema.validate(invalid_time)
    except SchemaError as e:
        assert "'Time' does not match 'Tijd'" in str(e)


@pytest.mark.parametrize(
    "name, file, expected_error",
    [
        (
            "CMIP",
            "CMIP6/snw_day_CanESM5_historical_r1i1p1f1_gn_19910101-20101231.nc",
            "'cf_time_dimension_schema' Missing key: 'units'",
        ),
        (
            "ECMWF",
            "ECMWF/era5_t2m-d2m_2024.nc",
            "'cf_time_dimension_schema' Missing keys: 'axis', 'standard_name', 'units'",
        ),
        (
            "NASA",
            "NASA/daymet_v4_tmax_annavg_hi_1987.nc",
            "Or('time', 'Time') did not validate '24-hour day based on local time'",
        ),
    ],
)
def test_invalid_real_data_dimensions(cassini, name, file, expected_error):
    ds = xr.open_dataset(cassini.fetch(file))
    if name == "ECMWF":
        ds = ds.rename_dims(
            {"latitude": "lat", "longitude": "lon", "valid_time": "time"}
        )

    dimensions = {
        "dimensions": {
            dim: ds[dim].attrs for dim in ds.dims if dim in ["time", "lat", "lon"]
        }
    }

    try:
        assert cf_dimensions_schema.validate(dimensions)
    except SchemaError as e:
        assert expected_error in str(e)


def test_valid_url():
    valid_urls = [
        "http://example.com",
        "https://example.com",
        "ftp://example.com",
        "http://localhost:8000",
        "https://www.example.com/path/to/resource",
    ]

    for url in valid_urls:
        assert url_validate(url) is not None, f"URL validation failed for: {url}"


def test_invalid_url():
    invalid_urls = [
        "htp://example.com",  # Typo in scheme
        "://example.com",  # Missing scheme
        "http:/example.com",  # Single slash after http
        "http://",  # Incomplete URL
        "example.com",  # No scheme
        "http://-example.com",  # Invalid domain start
    ]

    for url in invalid_urls:
        assert (
            url_validate(url) is None
        ), f"URL validation should have failed for: {url}"
