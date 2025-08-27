from __future__ import annotations

import pytest
import xarray as xr
from schema import SchemaError

from miranda.convert.corrections import CONFIG_FILES
from miranda.validate import (
    cf_dimensions_schema,
    cf_header_schema,
    cf_variables_schema,
    url_validate,
    validate_json,
)


class TestValidateJSON:
    @pytest.mark.xfail(raises=SchemaError, strict=False)
    @pytest.mark.parametrize("project, configuration", ([(k, v) for k, v in CONFIG_FILES.items()]))
    def test_converter_files(self, project, configuration):
        """Test that the converter configurations are valid JSON according to existing schemas."""
        try:
            validate_json(configuration)
        except ValueError as e:
            if project in [
                "CMIP",
                "CORDEX",
                "ESPO-G6-E5L",
                "ESPO-G6-R2",
                "NEX-GDDP-CMIP6",
            ]:
                # These projects have a different schema that is not yet implemented
                assert "Schema is not CF-compliant. No validation is possible." in str(e)
            else:
                raise e


class TestHeaderSchema:
    def test_valid_header_schema(self):
        valid_header_input = {
            "title": "Test Dataset",
            "institution": "Test Institution",
            "source": "Test Source",
            "history": "Created on 2023-10-01",
            "references": "https://example.com/references",
            "Conventions": "CF-1.8",
            "license": "CC-BY-4.0",
            "license_type": "proprietary",
            "processing_level": "raw",
            "table_id": "test-321",
            "type": "test_data",
            "_miranda_version": True,
            "_special_attrs": {
                "myproject": "1.8"
            },  # underscore keys that aren't "_freq" or "_miranda_version" must point to either a dictionary or a boolean
            "_project_name": {"myproject": True},
        }

        assert cf_header_schema.validate(valid_header_input)

    def test_invalid_header_schema(self):
        invalid_header_input = {
            "title": "Test Dataset",
            "institution": "Test Institution",
            "source": "Test Source",
            "history": "Created on 2023-10-01",
            "references": "https://example.com/references",
            "_miranda_version": True,
            "_special_attrs": {"myproject": "1.8"},
            "_project_name": "Test Project",  # underscore keys that aren't "_freq" or "_miranda_version" must point to a dictionary
        }
        try:
            cf_header_schema.validate(invalid_header_input)
        except SchemaError as e:
            assert "'Test Project' should be instance of 'dict'" in str(e)


class TestVariablesSchema:
    def test_valid_variables_schema(self):
        valid_variable_input = {
            "some-variable": {
                "standard_name": "air_temperature",
                "_cf_variable_name": "temperature",
                "_corrected_units": "K",
                "units": "K",
                "long_name": "Air Temperature",
            }
        }

        assert cf_variables_schema.validate(valid_variable_input)

    def test_invalid_variable_schema(self):
        invalid_variable_input = {
            "T": {
                "standard_name": "air_temperature",
                "_cf_variable_name": "temp",
                "_corrected_units": 300,  # Invalid type
                "units": "K",
            }
        }

        try:
            cf_variables_schema.validate(invalid_variable_input)
        except SchemaError as e:
            assert "300 should be instance of 'str'" in str(e)


class TestDimensionsSchema:
    def test_valid_dimensions_schema(self):
        valid_dimensions_input = {
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
        assert cf_dimensions_schema.validate(valid_dimensions_input)

    def test_invalid_time_schema(self):
        invalid_dimensions_time = {
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
        try:
            cf_dimensions_schema.validate(invalid_dimensions_time)
        except SchemaError as e:
            assert "'Tijd' does not match '^time$|^Time$'" in str(e)

    @pytest.mark.parametrize(
        "name, file, expected_error",
        [
            (
                "CMIP",
                "CMIP6/snw_day_CanESM5_historical_r1i1p1f1_gn_19910101-20101231.nc",
                "Time dimension must contain either 'units' or '_units', but not both",
            ),
            (
                "ECMWF",
                "ECMWF/era5_t2m-d2m_2024.nc",
                "'cf_time_dimension_schema' Missing keys: 'axis', 'standard_name'",
            ),
            (
                "NASA",
                "NASA/daymet_v4_tmax_annavg_hi_1987.nc",
                "'24-hour day based on local time' does not match '^time$|^Time$'",
            ),
        ],
    )
    def test_invalid_real_data_dimensions(self, cassini, name, file, expected_error):
        ds = xr.open_dataset(cassini.fetch(file))
        if name == "ECMWF":
            ds = ds.rename_dims({"latitude": "lat", "longitude": "lon", "valid_time": "time"})

        dimensions = {dim: ds[dim].attrs for dim in ds.dims if dim in ["time", "lat", "lon"]}

        try:
            assert cf_dimensions_schema.validate(dimensions)
        except SchemaError as e:
            assert expected_error in str(e)


class TestURLValidation:
    def test_valid_url(self):
        valid_urls = [
            "http://example.com",
            "https://example.com",
            "ftp://example.com",
            "http://localhost:8000",
            "https://www.example.com/path/to/resource",
        ]

        for url in valid_urls:
            assert url_validate(url) is not None, f"URL validation failed for: {url}"

    def test_invalid_url(self):
        invalid_urls = [
            "htp://example.com",  # Typo in scheme
            "://example.com",  # Missing scheme
            "http:/example.com",  # Single slash after http
            "http://",  # Incomplete URL
            "example.com",  # No scheme
            "http://-example.com",  # Invalid domain start
        ]

        for url in invalid_urls:
            assert url_validate(url) is None, f"URL validation should have failed for: {url}"
