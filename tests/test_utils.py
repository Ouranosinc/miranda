from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import pytest  # noqa

import miranda.eccc._utils as eccc_utils  # noqa
import miranda.utils


class TestWorkingDirectory:
    def test_different_directory(self):
        present_directory = Path.cwd()
        parent = present_directory.parent

        with miranda.utils.working_directory(parent) as location:
            assert location != present_directory
            assert Path.cwd() == parent
        assert Path.cwd() == present_directory


class TestEnvCanVariables:
    def test_hourly_cf_dictionaries(self):
        keys = [76, 77, 78, 80, 123, 262]

        codes = list()
        variables = dict()
        for key in keys:
            variables[key] = eccc_utils.cf_station_metadata(key)
            codes.append(variables[key]["standard_name"])
            if variables[key]["standard_name"] == "dry_bulb_temperature":
                assert variables[key]["raw_units"] == "degC"
                assert variables[key]["units"] == "K"
            assert variables[key]["missing_flags"] == "M"

        assert set(codes) == {
            "wind_speed_u2a",
            "atmospheric_pressure",
            "dry_bulb_temperature",
            "relative_humidity",
            "rainfall_flux",
            "precipitation_amount",
        }

    def test_daily_cf_dictionaries(self):
        keys = [
            1,
            2,
            3,
            10,
            11,
            12,
        ]

        codes = list()
        variables = dict()
        for key in keys:
            variables[key] = eccc_utils.cf_station_metadata(key)
            codes.append(variables[key]["standard_name"])
            if variables[key]["standard_name"].startswith("air_temperature"):
                assert variables[key]["raw_units"] == "degC"
                assert variables[key]["units"] == "K"
            elif variables[key]["standard_name"].endswith("precipitation_amount"):
                assert variables[key]["raw_units"] in ["cm", "mm"]
                assert variables[key]["units"] == "m"
            assert variables[key]["missing_flags"] == "M"

        assert set(codes) == {
            "air_temperature",
            "air_temperature_maximum",
            "air_temperature_minimum",
            "precipitation_amount",
            "liquid_precipitation_amount",
            "solid_precipitation_amount",
        }


class TestCreationDate:
    def test_newly_created_file(self):
        filename = "testfile.txt"
        testfile = Path.cwd().joinpath(filename)

        with open(testfile.name, "w") as f:
            f.write(filename)

        assert miranda.io.utils.creation_date(testfile) == date.today()
        testfile.unlink()


class TestReadPrivileges:
    def test_allowed_folder(self):
        here = Path.cwd()
        allowed = miranda.utils.read_privileges(here)
        assert allowed

    def test_nonexistent_folder_strict(self):
        mythical_folder = Path("/here/there/everywhere")
        with pytest.raises(OSError):
            miranda.utils.read_privileges(mythical_folder, strict=True)

    @pytest.mark.skipif(os.name != "posix", reason="not Windows")
    def test_forbidden_folder_lax(self):
        root_folder = Path.cwd().root
        if os.getenv("CI"):
            root_folder = Path(Path.cwd().root) / "root"
        allowed = miranda.utils.read_privileges(root_folder, strict=False)
        if os.getenv("CI"):
            assert not allowed
        else:
            assert allowed
