import logging
import multiprocessing as mp
import os
import re
import warnings
from functools import partial
from logging import config
from os import PathLike
from pathlib import Path
from types import GeneratorType
from typing import Dict, List, Optional, Union

import netCDF4 as nc  # noqa
import pandas as pd
import schema
import zarr
from pandas._libs.tslibs import NaTType  # noqa

from miranda.cv import INSTITUTIONS, PROJECT_MODELS
from miranda.scripting import LOGGING_CONFIG
from miranda.validators import FACETS_SCHEMA

from ._time import (
    TIME_UNITS_TO_FREQUENCY,
    TIME_UNITS_TO_TIMEDELTA,
    DecoderError,
    date_parser,
)

config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "Decoder",
    "guess_project",
]


def guess_project(file: Union[os.PathLike, str]) -> str:
    file_name = Path(file).stem

    potential_names = file_name.split("_")
    for project, models in PROJECT_MODELS.items():
        if any([model in potential_names for model in models]):
            return project
    raise DecoderError(f"Unable to determine project from file name: '{file_name}'.")


class Decoder:

    project = None
    _file_facets = dict()

    def __init__(self, project: Optional[str]):
        self.project = project

    @staticmethod
    def _decoder(
        d: dict,
        fail_early: bool,
        proj: str,
        lock,
        file: Union[str, Path],
    ) -> None:
        """

        Notes
        -----
        lock is a threading lock object
        """
        with lock:
            if proj is None:
                try:
                    proj = guess_project(file)
                except DecoderError:
                    print(
                        "Unable to determine 'activity': Signature for 'activity' must be set manually for file: "
                        f"{file}."
                    )
                    if fail_early:
                        raise

            decode_function_name = f"decode_{proj.lower().replace('-','_')}"
            try:
                _deciphered = getattr(Decoder, decode_function_name)(Path(file))
                if fail_early:
                    FACETS_SCHEMA.validate(_deciphered)
                print(
                    f"Deciphered the following from {Path(file).name}: {_deciphered.items()}"
                )
                d[file] = _deciphered
            except AttributeError as e:
                print(f"Unable to read data from {Path(file).name}: {e}")
            except schema.SchemaError as e:
                print(f"Decoded facets from {Path(file).name} are not valid: {e}")
                raise

    def decode(
        self,
        files: Union[os.PathLike, str, List[Union[str, os.PathLike]], GeneratorType],
        raise_error: bool = False,
    ) -> None:
        """Decode facets from file or list of files.

        Parameters
        ----------
        files: Union[str, Path, List[Union[str, Path]]]
        raise_error: bool
        """
        if isinstance(files, (str, os.PathLike)):
            files = [files]
        if self.project is None:
            warnings.warn(
                "The decoder 'project' is not set; Decoding step will be much slower."
            )
        else:
            logging.info(f"Deciphering metadata with project = '{self.project}'")
        manager = mp.Manager()
        _file_facets = manager.dict()
        lock = manager.Lock()
        func = partial(self._decoder, _file_facets, raise_error, self.project, lock)

        with mp.Pool() as pool:
            pool.imap(func, files, chunksize=10)
            pool.close()
            pool.join()

        self._file_facets.update(_file_facets)

    def facets_table(self):
        raise NotImplementedError()

    def file_facets(self) -> Dict[os.PathLike, Dict]:
        return self._file_facets

    @classmethod
    def _from_dataset(cls, file: Union[Path, str]) -> (str, str, Dict):
        file_name = Path(file).stem

        variable_name = cls._decode_primary_variable(file)
        variable_date = file_name.split("_")[-1]

        if file.is_file() and file.suffix in [".nc", ".nc4"]:
            ds = nc.Dataset(file)
            data = dict()
            for k in ds.ncattrs():
                data[k] = getattr(ds, k)
        elif file.is_dir() and file.suffix == ".zarr":
            ds = zarr.open(file, mode="r")
            data = ds.attrs.asdict()
        else:
            raise DecoderError(f"Unable to read dataset: `{file.name}`.")
        return variable_name, variable_date, data

    @staticmethod
    def _decode_primary_variable(file: Path) -> str:
        """Attempts to find the primary variable of a netCDF

        Parameters
        ----------
        file: Union[Path, str]

        Returns
        -------
        str
        """
        dimsvar_dict = dict()
        coords = ("time", "lat", "lon", "rlat", "rlon", "height", "lev", "rotated_pole")
        suggested_variable = file.name.split("_")[0]

        if file.is_file() and file.suffix in [".nc", ".nc4"]:
            data = nc.Dataset(file, mode="r")
            for var_name, var_attrs in data.variables.items():
                dimsvar_dict[var_name] = {
                    k: var_attrs.getncattr(k) for k in var_attrs.ncattrs()
                }
            for k in dimsvar_dict.keys():
                if not str(k).startswith(coords) and suggested_variable == k:
                    return str(k)

        elif file.is_dir() and file.suffix == ".zarr":
            data = zarr.open(str(file), mode="r")
            for k in data.array_keys():
                if not str(k).startswith(coords) and suggested_variable == k:
                    return str(k)
        else:
            raise NotImplementedError()

    @staticmethod
    def _decode_time_info(
        file: Optional[Union[PathLike, str, List[str]]] = None,
        data: Optional[Dict] = None,
        *,
        field: str,
    ) -> Union[str, NaTType]:
        """

        Parameters
        ----------
        file: Union[os.PathLike, str] or List[str], optional
        data: dict, optional
        field: {"timedelta", "frequency"}

        Returns
        -------

        """
        if not file and not data:
            raise ValueError()

        if field == "frequency":
            time_dictionary = TIME_UNITS_TO_FREQUENCY
        elif field == "timedelta":
            time_dictionary = TIME_UNITS_TO_TIMEDELTA
        else:
            raise NotImplementedError()

        if isinstance(file, (str, PathLike)):
            file = Path(file).name.split("_")

        if isinstance(file, list):
            potential_times = [segment in file for segment in time_dictionary.keys()]
            if potential_times:
                if potential_times[0] in ["fx", "fixed"]:
                    if field == "timedelta":
                        return pd.NaT
                    return "fx"
                if field == "timedelta":
                    return pd.to_timedelta(time_dictionary[potential_times[0]])
                return time_dictionary[potential_times[0]]
        elif data:
            potential_time = data["frequency"]
            if potential_time == "":
                time_units = data["time"].units
                potential_time = time_units.split()[0]
            if field == "timedelta":
                if potential_time in ["fx", "fixed"]:
                    return pd.NaT
                return pd.to_timedelta(time_dictionary[potential_time])
            return time_dictionary[potential_time]

    @classmethod
    def decode_converted(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        facets = dict()
        facets.update(data)
        del facets["history"]

        facets["date"] = date

        file_format = data.get("output_format")
        if format:
            facets["format"] = file_format
        else:
            facets["format"] = data["format"]

        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["variable"] = variable

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @staticmethod
    def decode_eccc_obs(self, file: Union[PathLike, str]) -> dict:
        raise NotImplementedError()

    @staticmethod
    def decode_ahccd_obs(self, file: Union[PathLike, str]) -> dict:
        raise NotImplementedError()

    @staticmethod
    def decode_melcc_obs(self, file: Union[PathLike, str]) -> dict:
        raise NotImplementedError()

    @classmethod
    def decode_pcic_candcs_u6(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        facets = dict()
        facets["activity"] = data["activity_id"]
        facets["mip_era"] = data["project_id"]
        facets["bias_adjust_institution"] = "PCIC"
        facets["date"] = date
        facets["domain"] = data["domain"]
        facets["experiment"] = str(data["GCM__experiment_id"]).replace(",", "-")
        facets["format"] = "netcdf"
        facets["frequency"] = data["frequency"]
        facets["institution"] = data["GCM__institution_id"]
        facets["member"] = (
            f"r{data['GCM__realization_index']}"
            f"i{data['GCM__initialization_index']}"
            f"p{data['GCM__physics_index']}"
            f"f{data['GCM__forcing_index']}"
        )
        facets["processing_level"] = "biasadjusted"
        facets["bias_adjust_project"] = "CanDCS-U6"
        facets["source"] = data["GCM__source_id"]
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable
        facets["version"] = data["GCM__data_specs_version"]

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_cmip6(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        facets = dict()
        facets["activity"] = data["activity_id"]
        facets["date"] = date
        facets["domain"] = "global"
        facets["experiment"] = data["experiment_id"]
        facets["format"] = "netcdf"
        facets["frequency"] = data["frequency"]
        facets["grid_label"] = data["grid_label"]
        facets["institution"] = data["institution_id"]
        facets["member"] = data["variant_label"]
        facets["modeling_realm"] = data["realm"]
        facets["processing_level"] = "raw"
        facets["mip_era"] = data["mip_era"]
        facets["source"] = data["source_id"]
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable

        try:
            facets["version"] = data["version"]
        except KeyError:
            possible_version = Path(file).parent.name
            if re.match(r"^[vV]\d+", possible_version):
                facets["version"] = Path(file).parent.name
            else:
                facets["version"] = "vNotFound"

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_cmip5(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        facets = dict()
        facets["activity"] = "CMIP"
        facets["date"] = date
        facets["domain"] = "global"
        facets["experiment"] = data["experiment_id"]
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(data=data, field="frequency")
        facets["institution"] = data["institute_id"]
        facets["member"] = data["parent_experiment_rip"]
        facets["modeling_realm"] = data["modeling_realm"]
        facets["processing_level"] = "raw"
        facets["mip_era"] = data["project_id"]
        facets["source"] = data["model_id"]
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable

        try:
            facets["version"] = data["version"]
        except KeyError:
            possible_version = Path(file).parent.name
            if re.match(r"^[vV]\d+", possible_version):
                facets["version"] = Path(file).parent.name
            else:
                facets["version"] = "vNotFound"

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_cordex(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        # FIXME: What to do about our internal data that breaks all established conventions?
        facets = dict()
        facets["activity"] = "CORDEX"

        if data.get("project_id") == "" or data.get("project_id") is None:
            facets["mip_era"] = "internal"
        elif data.get("project_id") == "CORDEX":
            facets["mip_era"] = "CMIP5"

        if date == "r0i0p0":
            facets["date"] = "fx"
        else:
            facets["date"] = date

        domain = data.get("CORDEX_domain").strip()
        if domain:
            facets["domain"] = domain
        else:
            domain = data.get("ouranos_domain_name").strip()
            if domain:
                facets["domain"] = domain
            else:
                msg = f"File {Path(file).name} has a nonstandard domain name."
                logging.error(msg)
                raise NotImplementedError(msg)

        # CORDEX-NAM on AWS mis-attributes the domain (22/44 should be 22i/44i)
        aws_keys = data.get("intake_esm_dataset_key")
        if aws_keys:
            facets["domain"] = aws_keys.split(".")[3]

        # The logic here is awful, but the information is bad to begin with.
        driving_institution_parts = str(data["driving_model_id"]).split("-")
        if driving_institution_parts[0] in INSTITUTIONS:
            driving_institution = driving_institution_parts[0]
        elif "-".join(driving_institution_parts[:2]) in INSTITUTIONS:
            driving_institution = "-".join(driving_institution_parts[:2])
        elif "-".join(driving_institution_parts[:3]) in INSTITUTIONS:
            driving_institution = "-".join(driving_institution_parts[:3])
        elif driving_institution_parts[0] == "GFDL":
            driving_institution = "NOAA-GFDL"
        elif driving_institution_parts[0] == "MPI":
            driving_institution = "MPI-M"
        else:
            raise AttributeError(
                "driving_institution (from driving_model_id: "
                f"`{data['driving_model_id']}`) is not valid."
            )

        facets["driving_institution"] = driving_institution
        if data["driving_model_id"].startswith("MPI"):
            facets["driving_model"] = f"MPI-M-{data['driving_model_id']}"
        else:
            facets["driving_model"] = data["driving_model_id"]
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(data=data, field="frequency")

        if data["institute_id"].strip() == "Our.":
            facets["institution"] = "Ouranos"
        else:
            facets["institution"] = data["institute_id"].strip()

        facets["processing_level"] = "raw"
        facets["source"] = data["model_id"]
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable

        try:
            facets["version"] = data["version"]
        except KeyError:
            possible_version = Path(file).parent.name
            if re.match(r"^[vV]\d+", possible_version):
                facets["version"] = Path(file).parent.name
            else:
                facets["version"] = "vNotFound"

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        try:
            facets["experiment"] = data["experiment_id"].strip()
        except KeyError:
            facets["experiment"] = data["driving_experiment_name"].strip()

        try:
            for potential_member in ["parent_experiment_rip", "parent_experiment"]:
                facets["member"] = data.get(potential_member)
                if facets["member"] == "N/A":
                    raise KeyError()
                else:
                    break
            if facets["member"] is None:
                raise KeyError()
        except KeyError:
            facets["member"] = data["driving_model_ensemble_member"].strip()

        return facets

    @classmethod
    def decode_isimip_ft(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        facets = dict()
        facets["activity"] = "ISIMIP"
        facets["mip_era"] = data["project_id"]

        facets["date"] = date
        facets["domain"] = "global"
        facets["co2_forcing_id"] = data["co2_forcing_id"]
        facets["experiment"] = data["experiment_id"]
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(data=data, field="frequency")
        facets["impact_model"] = data["impact_model_id"]
        facets["institution"] = data["institute_id"]
        facets["member"] = data["driving_model_ensemble_member"]
        facets["modeling_realm"] = data["modeling_realm"]
        facets["social_forcing_id"] = data["social_forcing_id"]
        facets["source"] = data["model_id"]
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets
