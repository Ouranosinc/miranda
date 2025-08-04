from __future__ import annotations

import pytest
import xarray as xr
from schema import SchemaError

from miranda.validate import cf_dimensions_schema


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
            },
            "lat": {
                "_cf_dimension_name": "lat",
                "_precision": 2,
                "axis": "Y",
                "standard_name": "latitude",
                "long_name": "Latitude",
            },
            "lon": {
                "_cf_dimension_name": "lon",
                "_precision": 2,
                "axis": "X",
                "standard_name": "longitude",
                "long_name": "Longitude",
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
