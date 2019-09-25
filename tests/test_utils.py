from datetime import date
from pathlib import Path

import pytest

from miranda import utils


class TestWorkingDirectory:
    def test_different_directory(self):
        present_directory = Path.cwd()
        parent = present_directory.parent

        with utils.working_directory(parent) as location:
            assert location != present_directory
            assert Path.cwd() == parent
        assert Path.cwd() == present_directory


class TestEnvCanVariables:
    def test_variable_dictionaries(self):
        keys = [
            "wind_speed",
            "station_pressure",
            "dry_bulb_temperature",
            "relative_humidity",
            "hourly_rainfall",
            "precipitation_amount",
        ]

        codes = list()
        variables = dict()
        for key in keys:
            variables[key] = utils.get_info_var(key)
            codes.append(variables[key]["code_var"])
            assert variables[key]["fact_add"] == 0
            assert variables[key]["flag_manquants"] == ["M"]
            assert variables[key]["least_significant_digit"] is None

        assert codes == [76, 77, 78, 80, 123, 262]


class TestCreationDate:
    def test_newly_created_file(self):
        filename = "testfile.txt"
        testfile = Path.cwd().joinpath(filename)

        with open(testfile.name, "w") as f:
            f.write(filename)

        assert utils.creation_date(testfile) == date.today()
        testfile.unlink()


class TestFolderOperations:
    def test_local_folder_creation(self):
        foldername = "testfolder"
        testfolder = Path.cwd().joinpath(foldername)

        utils.make_local_dirs(testfolder)
        assert testfolder.exists()
        assert testfolder.stat().st_mode == 0o40775
        testfolder.rmdir()

    def test_make_forbidden_folder(self):
        foldername = "testfolder"
        root_drive = Path.cwd().root
        testfolder = Path(root_drive).joinpath(foldername)

        with pytest.raises(OSError):
            utils.make_local_dirs(testfolder)
