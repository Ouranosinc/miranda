from __future__ import annotations

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

import netCDF4 as nc  # noqa
import pandas as pd
import schema
import xarray as xr
import zarr
from pandas._libs.tslibs import NaTType  # noqa

from miranda.convert.utils import date_parser, find_version_hash  # noqa
from miranda.cv import VALIDATION_ENABLED
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

from ._time import TIME_UNITS_TO_FREQUENCY, TIME_UNITS_TO_TIMEDELTA, DecoderError

if VALIDATION_ENABLED:
    from miranda.cv import INSTITUTIONS, PROJECT_MODELS
    from miranda.validators import FACETS_SCHEMA  # noqa


config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "Decoder",
    "guess_project",
]


def guess_project(file: os.PathLike | str) -> str:
    """Guess the name of the project

    Parameters
    ----------
    file : str or os.PathLike

    Returns
    -------
    str
    """
    file_name = Path(file).stem

    potential_names = file_name.split("_")
    if VALIDATION_ENABLED:
        for project, models in PROJECT_MODELS.items():
            if any([model in potential_names for model in models]):
                return project
        raise DecoderError(
            f"Unable to determine project from file name: '{file_name}'."
        )
    raise DecoderError("Project determination requires pyessv-archive source files.")


class Decoder:
    project = None
    guess = False
    _file_facets = dict()

    def __init__(self, project: str | None):
        self.project = project

    @staticmethod
    def _decoder(
        d: dict,
        fail_early: bool,
        proj: str,
        guess: bool,
        lock: mp.Lock,
        file: str | Path,
    ) -> None:
        if proj is None:
            if guess:
                try:
                    proj = guess_project(file)
                except DecoderError:
                    print(
                        "Unable to determine 'activity': Signature for 'activity' must be set manually for file: "
                        f"{file}."
                    )
                    if fail_early:
                        raise
            else:
                proj = "converted"

        decode_function_name = f"decode_{proj.lower().replace('-','_')}"
        try:
            with lock:
                _deciphered = getattr(Decoder, decode_function_name)(Path(file))
                if fail_early:
                    if VALIDATION_ENABLED:
                        FACETS_SCHEMA.validate(_deciphered)
                    else:
                        print(
                            "Validation requires pyessv-archive source files. Skipping validation checks."
                        )
                print(
                    f"Deciphered the following from {Path(file).name}:\n"
                    f"{_deciphered.items()}"
                )
                d[file] = _deciphered

        except (AttributeError, NotImplementedError):
            print(f"Unable to read data from {Path(file)}. Ensure pathname is correct.")
            raise
        except schema.SchemaError as e:
            print(f"Decoded facets from {Path(file).name} are not valid: {e}")

    def decode(
        self,
        files: os.PathLike | str | list[str | os.PathLike] | GeneratorType,
        chunks: int | None = None,
        raise_error: bool = False,
    ) -> None:
        """Decode facets from file or list of files.

        Parameters
        ----------
        files : str or Path or list of str or Path or generator
        chunks : int, optional
            The chunk size used when processing files. Not to be confused with xarray chunks for dimensions.
        raise_error : bool
        """
        if isinstance(files, (str, os.PathLike)):
            files = [files]

        if chunks is None and isinstance(files, list):
            if len(files) >= 10:
                chunk_size = 10
            elif 1 <= len(files) < 10:
                chunk_size = len(files)
            else:
                raise ValueError("No file entries found.")
        elif isinstance(files, GeneratorType):
            chunk_size = 10
        else:
            chunk_size = chunks

        if self.project is None:
            warnings.warn(
                "The decoder 'project' is not set; Decoding step will be much slower."
            )
        else:
            logging.info(f"Deciphering metadata with project = '{self.project}'")

        with mp.Manager() as manager:
            _file_facets = manager.dict()
            lock = manager.Lock()
            func = partial(
                self._decoder, _file_facets, raise_error, self.project, self.guess, lock
            )

            with mp.Pool() as pool:
                pool.imap(func, files, chunksize=chunk_size)
                pool.close()
                pool.join()

            self._file_facets.update(_file_facets)

    def facets_table(self):
        raise NotImplementedError()

    def file_facets(self) -> dict[os.PathLike, dict]:
        return self._file_facets

    @classmethod
    def _from_dataset(cls, file: Path | str) -> (str, str, dict):
        file_name = Path(file).stem

        try:
            variable_name = cls._decode_primary_variable(file)
        except DecoderError:
            logging.error(f"Unable to open dataset: {file.name}")
            raise

        datetimes = file_name.split("_")[-1]

        if file.is_file() and file.suffix in [".nc", ".nc4"]:
            with nc.Dataset(file, mode="r") as ds:
                data = dict()
                for k in ds.ncattrs():
                    data[k] = getattr(ds, k)
        elif file.is_dir() and file.suffix == ".zarr":
            with zarr.open(file, mode="r") as ds:
                data = ds.attrs.asdict()
        else:
            raise DecoderError(f"Unable to read dataset: `{file.name}`.")
        return variable_name, datetimes, data

    @staticmethod
    def _decode_primary_variable(file: Path) -> str:
        """Attempts to find the primary variable of a netCDF

        Parameters
        ----------
        file: Path

        Returns
        -------
        str
        """
        dimsvar_dict = dict()
        coords = (
            "height",
            "lat",
            "latitude",
            "lev",
            "level",
            "lon",
            "longitude",
            "rlat",
            "rlon",
            "rotated_pole",
            "time",
        )
        try:
            if file.is_file() and file.suffix in [".nc", ".nc4"]:
                with nc.Dataset(file, mode="r") as ds:
                    for var_name, var_attrs in ds.variables.items():
                        dimsvar_dict[var_name] = {
                            k: var_attrs.getncattr(k) for k in var_attrs.ncattrs()
                        }
                for k in dimsvar_dict.keys():
                    if not str(k).startswith(coords) and k in file.stem:
                        return str(k)

            elif file.is_dir() and file.suffix == ".zarr":
                with zarr.open(str(file), mode="r") as ds:
                    for k in ds.array_keys():
                        if not str(k).startswith(coords) and k in file.stem:
                            return str(k)
            else:
                raise NotImplementedError()
        except ValueError:
            raise DecoderError()

    @staticmethod
    def _decode_hour_of_day_info(
        file: PathLike | str,
    ) -> dict:
        """

        Parameters
        ----------
        file : Path or str

        Returns
        -------
        dict
        """
        if isinstance(file, str):
            file = Path(file)

        if file.is_file() and file.suffix in [".nc", ".nc4"]:
            with nc.Dataset(file, mode="r") as ds:
                if "time" in ds.variables.keys():
                    hour = nc.num2date(
                        ds["time"][0], ds["time"].units, ds["time"].calendar
                    ).hour
                else:
                    hour = None
            return dict(hour_of_day=hour)

        elif file.is_dir() and file.suffix == ".zarr":
            warnings.warn("This is not currently implemented")

            # with zarr.open(str(file), mode="r") as ds:
            #     if "time" in ds.array_keys():
            #         pass

            return dict()

        else:
            raise NotImplementedError()

    @staticmethod
    def _decode_time_info(
        file: PathLike | str | list[str] | None = None,
        data: dict | None = None,
        term: str | None = None,
        *,
        field: str = None,
    ) -> str | NaTType:
        """

        Parameters
        ----------
        file : os.PathLike or str, optional
        data : dict, optional
        term : str
        field : {"timedelta", "frequency"}

        Returns
        -------
        str or NaTType
        """
        if not file and not data and not term:
            raise ValueError("Nothing passed to parse time info from.")

        if field == "frequency":
            time_dictionary = TIME_UNITS_TO_FREQUENCY
        elif field == "timedelta":
            time_dictionary = TIME_UNITS_TO_TIMEDELTA
        else:
            raise NotImplementedError()

        if term:
            if term in ["fx", "fixed"]:
                if field == "timedelta":
                    return pd.NaT
                return "fx"
            return pd.to_timedelta(time_dictionary[term])

        if data and not file:
            potential_time = data.get("frequency")
            if not potential_time:
                if hasattr(data, "time"):
                    time_units = data["time"].units
                    potential_time = time_units.split()[0]
                else:
                    logging.warning(
                        f"Could not find `frequency` or `time` for {Path(file).name}. Assuming `fx`."
                    )
                    potential_time = "fx"
            if potential_time in ["ymon", "yseas", "fixed", "fx"]:
                logging.warning(f"Found `{potential_time}`. Frequency is likely `fx`.")
                if field == "frequency":
                    return "fx"
                if field == "timedelta":
                    return pd.NaT
                raise ValueError()

            if field == "timedelta":
                if potential_time in ["fx", "fixed"]:
                    return pd.NaT
                return pd.to_timedelta(time_dictionary[potential_time])
            return time_dictionary[potential_time]

        if file and not data:
            for delimiter in ["_", "."]:
                file_parts = Path(file).stem.split(delimiter)
                potential_times = [
                    segment
                    for segment in file_parts
                    if segment in time_dictionary.keys()
                ]
                if potential_times:
                    if potential_times[0] in ["fx", "fixed"]:
                        if field == "frequency":
                            return "fx"
                        if field == "timedelta":
                            return pd.NaT
                        raise ValueError(f"Field `{field}` not supported.")
                    if field == "timedelta":
                        return pd.to_timedelta(time_dictionary[potential_times[0]])
                    return time_dictionary[potential_times[0]]

        if file and data:
            for delimiter in ["_", "."]:
                file_parts = Path(file).stem.split(delimiter)
                potential_times = [
                    segment
                    for segment in file_parts
                    if segment in time_dictionary.keys()
                ]
                potential_time = data.get("frequency", "")
                if potential_time == "":
                    if hasattr(data, "time"):
                        time_units = data["time"].units
                        potential_time = time_units.split()[0]
                    else:
                        logging.warning(
                            f"Could not find `frequency` or `time` for {Path(file).name}. Assuming `fx`."
                        )
                        potential_time = "fx"
                if potential_time in ["ymon", "yseas", "fixed", "fx"]:
                    logging.warning(
                        f"Found `{potential_time}`. Frequency is likely `fx`."
                    )
                    if "fx" in file_parts or "fixed" in file_parts:
                        if field == "frequency":
                            return "fx"
                        if field == "timedelta":
                            return pd.NaT
                        raise ValueError(f"Field `{field}` not supported.")

                if potential_time in potential_times:
                    return time_dictionary[potential_time]
                elif potential_times:
                    break

            logging.warning(
                f"Frequency from metadata (`{potential_time}`) not found in filename (`{Path(file).name}`): "
                "Performing more rigorous frequency checks."
            )
            if Path(file).is_file() and Path(file).suffix in [".nc", ".nc4"]:
                engine = "netcdf4"
            elif Path(file).is_dir() and Path(file).suffix == ".zarr":
                engine = "zarr"
            else:
                raise DecoderError(
                    f"File is not valid netcdf or zarr: {Path(file).name}"
                )

            _ds = xr.open_dataset(
                file,
                engine=engine,
                drop_variables="time_bnds",
            )
            if not hasattr(_ds, "time"):
                logging.warning(
                    "Dataset does not contain time array. Assuming fixed variable."
                )
                if field == "frequency":
                    return "fx"
                if field == "timedelta":
                    return pd.NaT
                raise ValueError(f"Field `{field}` not supported.")
            else:
                _, found_freq = get_time_frequency(_ds.time)

            if found_freq in potential_times:
                logging.warning(
                    "Time frequency found in dataset on analysis was found in filename. "
                    f"Metadata for `{Path(file).name} is probably incorrect. "
                    f"Basing fields on `{found_freq}`."
                )
                return time_dictionary[found_freq]
            elif found_freq in ["month", "mon"]:
                for f in ["Amon", "Omon", "monC", "monthly", "months", "mon"]:
                    if f in potential_times:
                        logging.warning(
                            "Month-like time frequency found in dataset on analysis was found in filename. "
                            f"Basing fields on `{f}`."
                        )
                        return time_dictionary[f]
            else:
                logging.warning(
                    "Time frequency found in dataset on analysis was not found in filename. "
                    f"Basing fields on `{found_freq}`."
                )
                return time_dictionary[found_freq]
        raise DecoderError(f"Time frequency indiscernible for file `{file}`.")

    @staticmethod
    def _decode_version(file: PathLike | str, data: dict) -> dict:
        """

        Parameters
        ----------
        file: os.PathLike or str
        data: dict

        Returns
        -------
        dict
        """
        version_info = dict()
        try:
            version_info["version"] = data["version"]
        except KeyError:
            possible_version = Path(file).parent
            if re.match(r"^[vV]\d+", possible_version.name):
                version_info["version"] = possible_version.name
            else:
                possible_version_signature = possible_version.glob(
                    f"{Path(file).stem}.v*"
                )
                for sig in possible_version_signature:
                    found_version = re.match(r"([vV]\d+)$", sig.suffix)
                    if found_version:
                        version_info["version"] = found_version.group()
                        version_info["sha256sum"] = sig.open().read()
                        break
                else:
                    version_info["version"] = "vNotFound"
        return version_info

    @classmethod
    def decode_converted(cls, file: PathLike | str) -> dict:
        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return facets

        facets.update(data)
        del facets["history"]

        facets["date"] = date

        file_format = data.get("output_format")
        if file_format:
            facets["format"] = file_format
        elif "format" in data:
            facets["format"] = data["format"]
        elif Path(file).suffix in [".nc", ".nc4"]:
            facets["format"] = "nc"
        elif Path(file).suffix in [".zarr"]:
            facets["format"] = "zarr"
        facets["variable"] = variable

        facets.update(cls._decode_version(data=data, file=file))
        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            if "frequency" not in facets:
                facets["timedelta"] = cls._decode_time_info(
                    data=data, file=file, field="frequency"
                )
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @staticmethod
    def decode_eccc_obs(self, file: PathLike | str) -> dict:
        raise NotImplementedError()

    @staticmethod
    def decode_ahccd_obs(self, file: PathLike | str) -> dict:
        raise NotImplementedError()

    @staticmethod
    def decode_melcc_obs(self, file: PathLike | str) -> dict:
        raise NotImplementedError()

    @classmethod
    def decode_pcic_candcs_u6(cls, file: PathLike | str) -> dict:
        if "Derived" in Path(file).parents:
            raise NotImplementedError("Derived CanDCS-U6 variables are not supported.")

        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return facets

        facets["activity"] = data["activity_id"]
        facets["mip_era"] = data["project_id"]
        facets["bias_adjust_institution"] = "PCIC"
        facets["date"] = date
        facets["domain"] = data["domain"]
        facets["experiment"] = str(data["GCM__experiment_id"]).replace(",", "-")
        facets["format"] = "netcdf"
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
        facets["type"] = "simulation"
        facets["variable"] = variable

        facets["version"] = f"v{data.get('GCM__data_specs_version')}"
        if facets["version"] is None:
            facets.update(find_version_hash(file=file))

        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            facets["frequency"] = cls._decode_time_info(
                data=data, file=file, field="frequency"
            )
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_cmip6(cls, file: PathLike | str) -> dict:
        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return facets

        facets["activity"] = data["activity_id"]
        facets["date"] = date
        facets["domain"] = "global"
        facets["experiment"] = data["experiment_id"]
        facets["format"] = "netcdf"
        facets["grid_label"] = data["grid_label"]
        facets["institution"] = data["institution_id"]
        facets["member"] = data["variant_label"]
        facets["modeling_realm"] = data["realm"]
        facets["processing_level"] = "raw"
        facets["mip_era"] = data["mip_era"]
        facets["source"] = data["source_id"]
        facets["type"] = "simulation"
        facets["variable"] = variable
        facets.update(cls._decode_version(data=data, file=file))
        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            facets["frequency"] = cls._decode_time_info(
                data=data, file=file, field="frequency"
            )
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_cmip5(cls, file: PathLike | str) -> dict:
        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return facets

        facets["activity"] = "CMIP"
        facets["date"] = date
        facets["domain"] = "global"
        facets["experiment"] = data["experiment_id"]
        facets["format"] = "netcdf"
        facets["institution"] = data["institute_id"]
        facets["member"] = data["parent_experiment_rip"]
        facets["modeling_realm"] = data["modeling_realm"]
        facets["processing_level"] = "raw"
        facets["mip_era"] = data["project_id"]
        facets["source"] = data["model_id"]
        facets["type"] = "simulation"
        facets["variable"] = variable
        facets.update(cls._decode_version(data=data, file=file))
        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            facets["frequency"] = cls._decode_time_info(
                data=data, file=file, field="frequency"
            )
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_cordex(cls, file: PathLike | str) -> dict:
        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return dict()

        # FIXME: What to do about our internal data that breaks all established conventions?
        facets["activity"] = "CORDEX"

        if data.get("project_id") == "" or data.get("project_id") is None:
            facets["mip_era"] = "internal"
        elif data.get("project_id") == "CORDEX":
            facets["mip_era"] = "CMIP5"

        if date == "r0i0p0":
            facets["date"] = "fx"
        else:
            facets["date"] = date

        domain = data.get("CORDEX_domain")
        if domain:
            facets["domain"] = domain.strip()
        else:
            domain = data.get("ouranos_domain_name")
            if domain:
                facets["domain"] = domain.strip()
            else:
                msg = f"File {Path(file).name} has a nonstandard domain name."
                logging.error(msg)
                raise NotImplementedError(msg)

        # CORDEX-NAM on AWS mis-attributes the domain (22/44 should be 22i/44i)
        aws_keys = data.get("intake_esm_dataset_key")
        if aws_keys:
            facets["domain"] = aws_keys.split(".")[3]

        title = data.get("title")
        if title:
            regridded_domain_found = re.search(r"\w{3}-\d{2}i", title)
            if regridded_domain_found:
                facets["domain"] = regridded_domain_found.group()

        # The logic here is awful, but the information is bad to begin with.
        driving_model = ""
        driving_institution = ""

        driving_institution_parts = str(data["driving_model_id"]).split("-")
        if VALIDATION_ENABLED:
            if driving_institution_parts[0] in INSTITUTIONS:
                driving_institution = driving_institution_parts[0]
            elif "-".join(driving_institution_parts[:2]) in INSTITUTIONS:
                driving_institution = "-".join(driving_institution_parts[:2])
            elif "-".join(driving_institution_parts[:3]) in INSTITUTIONS:
                driving_institution = "-".join(driving_institution_parts[:3])
        else:
            logging.warning(
                "CORDEX Metadata validation checks require PyESSV. "
                "Driving institution cannot be determined."
            )
            driving_model = data["driving_model_id"]

        if data["driving_model_id"].startswith("GFDL"):
            driving_institution = "NOAA-GFDL"
            driving_model = f"NOAA-GFDL-{data['driving_model_id']}"
        elif data["driving_model_id"].startswith("MPI-ESM"):
            driving_institution = "MPI-M"
            driving_model = f"MPI-M-{data['driving_model_id']}"
        elif data["driving_model_id"].startswith("HadGEM2"):
            driving_institution = "MOHC"
            driving_model = f"MOHC-{data['driving_model_id']}"
        elif data["driving_model_id"].startswith("CNRM-CM5"):
            driving_institution = "CNRM-CERFACS"
            driving_model = f"CNRM-CERFACS-{data['driving_model_id']}"

        elif VALIDATION_ENABLED and not driving_institution:
            raise DecoderError(
                "driving_institution (from driving_model_id: "
                f"`{data['driving_model_id']}`) is not valid."
            )

        facets["driving_institution"] = driving_institution.strip()
        if driving_model:
            facets["driving_model"] = driving_model.strip()
        else:
            facets["driving_model"] = str(data["driving_model_id"]).strip()

        facets["format"] = "netcdf"

        if data["institute_id"].strip() == "Our.":
            facets["institution"] = "Ouranos"
        else:
            facets["institution"] = data["institute_id"].strip()

        facets["processing_level"] = "raw"
        facets["source"] = data["model_id"]
        facets["type"] = "simulation"
        facets["variable"] = variable

        facets.update(cls._decode_version(data=data, file=file))
        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            facets["frequency"] = cls._decode_time_info(
                data=data, file=file, field="frequency"
            )
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
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
    def decode_isimip_ft(cls, file: PathLike | str) -> dict:
        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return facets

        facets["activity"] = "ISIMIP"
        facets["mip_era"] = data["project_id"]
        facets["date"] = date
        facets["domain"] = "global"
        facets["co2_forcing_id"] = data["co2_forcing_id"]
        facets["experiment"] = data["experiment_id"]
        facets["format"] = "netcdf"
        facets["impact_model"] = data["impact_model_id"]
        facets["institution"] = data["institute_id"]
        facets["member"] = data["driving_model_ensemble_member"]
        facets["modeling_realm"] = data["modeling_realm"]
        facets["social_forcing_id"] = data["social_forcing_id"]
        facets["source"] = data["model_id"]
        facets["type"] = "simulation"
        facets["variable"] = variable

        facets.update(cls._decode_version(data=data, file=file))
        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            facets["frequency"] = cls._decode_time_info(data=data, field="frequency")
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_nex_gddp_cmip6(cls, file: PathLike | str) -> dict:
        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return facets

        facets["experiment"] = data["scenario"]
        facets["activity"] = (
            "CMIP" if facets["experiment"] == "historical" else "ScenarioMIP"
        )
        facets["institution"] = data["cmip6_institution_id"]
        facets["member"] = data["variant_label"]
        facets["processing_level"] = "biasadjusted"
        facets["bias_adjust_project"] = "NEX-GDDP-CMIP6"
        facets["bias_adjust_institution"] = "NASA"
        facets["mip_era"] = "CMIP6"
        facets["source"] = data["cmip6_source_id"]
        facets["type"] = "simulation"
        facets["variable"] = variable
        facets.update(cls._decode_version(data=data, file=file))
        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            facets["frequency"] = cls._decode_time_info(
                data=data, file=file, field="frequency"
            )
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_espo_g6_r2(cls, file: PathLike | str) -> dict:
        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return facets

        facets["bias_adjust_project"] = "ESPO-G6-R2"
        facets["processing_level"] = "biasadjusted"
        facets["version"] = "1.0.0"
        facets["domain"] = "NAM"
        for f in [
            "experiment",
            "activity",
            "institution",
            "member",
            "bias_adjust_institution",
            "mip_era",
            "source",
            "type",
        ]:
            facets[f] = data[f"cat:{f}"]
        facets["variable"] = variable
        # facets.update(cls._decode_version(data=data, file=file))
        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            facets["frequency"] = cls._decode_time_info(
                data=data, file=file, field="frequency"
            )
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets

    @classmethod
    def decode_espo_g6_e5l(cls, file: PathLike | str) -> dict:
        facets = dict()
        try:
            variable, date, data = cls._from_dataset(file=file)
        except DecoderError:
            return facets

        facets["bias_adjust_project"] = "ESPO-G6-E5L"
        facets["processing_level"] = "biasadjusted"
        facets["version"] = "1.0.0"
        facets["domain"] = "NAM"
        for f in [
            "experiment",
            "activity",
            "institution",
            "member",
            "bias_adjust_institution",
            "mip_era",
            "source",
            "type",
        ]:
            facets[f] = data[f"cat:{f}"]
        facets["variable"] = variable
        # facets.update(cls._decode_version(data=data, file=file))
        facets.update(cls._decode_hour_of_day_info(file=file))

        try:
            facets["frequency"] = cls._decode_time_info(
                data=data, file=file, field="frequency"
            )
            facets["timedelta"] = cls._decode_time_info(
                term=facets["frequency"], field="timedelta"
            )
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderError:
            pass

        return facets
