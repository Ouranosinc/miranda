import logging
import os
import re
import warnings
from logging import config
from os import PathLike
from pathlib import Path
from typing import Dict, List, Optional, Union

import netCDF4 as nc  # noqa
import pandas as pd
import schema
import zarr
from pandas._libs.tslibs import NaTType  # noqa

from miranda.scripting import LOGGING_CONFIG

from ._models import CMIP5_GCM_PROVIDERS, CMIP6_GCM_PROVIDERS, PROJECT_MODELS
from ._utils import DecoderException, date_parser

config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "Decoder",
    "guess_project",
]

BASIC_DT_VALIDATION = r"\s*(?=\d{2}(?:\d{2})?)"
DATE_VALIDATION = r"^\d{4}-(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])$"

TIME_UNITS_TO_FREQUENCY = {
    "subhrPt": "sub-hr",
    "hourly": "1hr",
    "hours": "1hr",
    "hr": "1hr",
    "6-hourly": "6hr",
    "daily": "day",
    "days": "day",
    "day": "day",
    "weekly": "sem",
    "weeks": "sem",
    "sem": "sem",
    "monthly": "mon",
    "months": "mon",
    "mon": "mon",
    "monC": "monC",
    "yearly": "yr",
    "years": "yr",
    "annual": "yr",
    "yr": "yr",
    "yrPt": "yrPt",
    "decadal": "dec",
    "decades": "dec",
    "dec": "dec",
    "fixed": "fx",
}

TIME_UNITS_TO_TIMEDELTA = {
    "hourly": "1h",
    "hours": "1h",
    "6-hourly": "6h",
    "daily": "1d",
    "day": "1d",
    "days": "1d",
    "weekly": "7d",
    "weeks": "7d",
    "sem": "7d",
    "mon": "30d",
    "monC": "30d",
    "QS": "90d",
    "qtr": "90d",
    "yearly": "365d",
    "years": "365d",
    "year": "365d",
}

facet_schema = schema.Schema(
    {
        schema.Optional("project"): str,
        "activity": str,
        "institution": str,
        "source": str,
        schema.Optional("driving_institution"): str,
        schema.Optional("driving_model"): str,
        schema.Optional("experiment"): str,
        "frequency": schema.And(
            str,
            lambda f: f
            in ["1hr", "3hr", "6hr", "day", "sem", "mon", "yr", "dec", "fx"],
        ),
        "domain": str,
        schema.Optional("member"): str,
        schema.Optional("variable"): str,
        schema.Optional("timedelta"): schema.Or(pd.Timedelta, NaTType),
        schema.Optional("date"): schema.Or(
            schema.Regex(BASIC_DT_VALIDATION, flags=re.I), "fx"
        ),
        schema.Optional("date_start"): schema.Or(
            schema.Regex(DATE_VALIDATION, flags=re.I), NaTType
        ),
        schema.Optional("date_end"): schema.Or(
            schema.Regex(DATE_VALIDATION, flags=re.I), NaTType
        ),
        schema.Optional("processing_level"): schema.And(
            str, lambda f: f in ["raw", "biasadjusted"]
        ),
        "format": schema.And(str, lambda f: f in ["netcdf", "zarr"]),
        schema.Optional("version"): str,
        "type": schema.And(
            str,
            lambda f: f
            in ["simulation", "reanalysis", "forecast", "gridded-obs", "station-obs"],
        ),
    },
    ignore_extra_keys=True,
)


def guess_project(file: Union[Path, str]) -> str:
    file_name = Path(file).stem

    potential_names = file_name.split("_")
    for project, models in PROJECT_MODELS.items():
        if any([model in potential_names for model in models]):
            return project
    raise DecoderException(f"Unable to determine project from file name: {file_name}.")


class Decoder:

    project = None
    _file_facets = dict()

    def __init__(self, project: Optional[str]):
        self.project = project

    def decode(
        self,
        files: Union[os.PathLike, str, List[Union[str, os.PathLike]]],
        raise_error=False,
    ):
        if isinstance(files, (str, os.PathLike)):
            files = [files]
        if self.project is None:
            warnings.warn(
                "The decoder 'project' is not set; Decoding step will be much slower."
            )

        _file_facets = dict()
        for file in files:
            if self.project is None:
                try:
                    project = guess_project(file)
                except DecoderException:
                    logging.error(
                        f"Signature for 'project' must be set manually for file: {file}."
                    )
                    if raise_error:
                        raise
                    else:
                        continue
            else:
                project = self.project

            _file_facets[file] = getattr(self, f"decode_{project.lower()}_netcdf")(
                Path(file)
            )
        self._file_facets.update(_file_facets)

    def facets_table(self):
        raise NotImplementedError()

    def file_facets(self) -> Dict[os.PathLike, Dict]:
        return self._file_facets

    @staticmethod
    def guess_project(file: Union[Path, str]) -> str:
        file_name = Path(file).stem

        potential_names = file_name.split("_")
        for project, models in PROJECT_MODELS.items():
            if any([model in potential_names for model in models]):
                return project
        raise DecoderException("Unable to determine project from file name.")

    @classmethod
    def _from_dataset(cls, file: Union[Path, str]) -> (str, str, nc.Dataset):
        file_name = Path(file).stem

        variable_name = cls._decode_primary_variable(file)
        variable_date = file_name.split("_")[-1]
        data = nc.Dataset(file)
        return variable_name, variable_date, data

    @staticmethod
    def _from_filename(file: Union[Path, str]) -> List[str]:
        file_name = Path(file).stem
        decode_file = file_name.split("_")
        return decode_file

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
        data: Optional[nc.Dataset] = None,
        *,
        field: str,
    ) -> Union[str, NaTType]:
        """

        Parameters
        ----------
        file
        data
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
            frequency = data.frequency
            if frequency == "":
                time_units = data["time"].units
                frequency = time_units.split()[0]
            if frequency in ["fx", "fixed", "mon", "Amon", "Emon"]:
                if field == "timedelta":
                    return pd.NaT
                return "fx"
            if field == "timedelta":
                return pd.to_timedelta(time_dictionary[frequency])
            if frequency in ["day"]:
                return frequency
            return time_dictionary[frequency]

    @staticmethod
    def decode_era5(file: Union[PathLike, str]) -> dict:
        raise NotImplementedError()

    @staticmethod
    def decode_generic_reanalysis(self, file: Union[PathLike, str]) -> dict:
        raise NotImplementedError()

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
    def decode_cmip6_netcdf(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        facets = dict()
        facets["activity"] = data.activity_id
        facets["date"] = date
        facets["domain"] = "global"
        facets["experiment"] = data.experiment_id
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(data=data, field="frequency")
        facets["institution"] = data.institution_id
        facets["member"] = data.variant_label
        facets["modeling_realm"] = data.realm
        facets["processing_level"] = "raw"
        facets["project"] = data.project
        facets["source"] = data.source_id
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable
        facets["version"] = data.version

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderException:
            pass

        facet_schema.validate(facets)

        logging.info(f"Deciphered the following from {file}: {facets.items()}")

        return facets

    @classmethod
    def decode_cmip6_name(cls, file: Union[PathLike, str]) -> dict:
        decode_file = cls._from_filename(file=file)

        facets = dict()
        facets["activity"] = "CMIP6"
        facets["date"] = decode_file[-1]
        facets["domain"] = "global"
        facets["experiment"] = decode_file[3]
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(file=decode_file, field="frequency")
        facets["grid_label"] = decode_file[5]
        facets["member"] = decode_file[4]
        facets["processing_level"] = "raw"
        facets["project"] = "CMIP6"
        facets["source"] = decode_file[2]
        facets["timedelta"] = cls._decode_time_info(file=decode_file, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = decode_file[0]

        if "mon" in facets["frequency"]:
            facets["frequency"] = "mon"

        try:
            logging.warning(
                "Using model to guess institute. Results may not be accurate."
            )
            facets["institution"] = CMIP6_GCM_PROVIDERS[facets["source"]]
        except KeyError:
            logging.info(f"Unable to find Institute for model: {facets['source']}")

        try:
            facets["date_start"] = date_parser(decode_file[-1])
            facets["date_end"] = date_parser(decode_file[-1], end_of_period=True)
        except DecoderException:
            pass

        facet_schema.validate(facets)

        logging.info(f"Deciphered the following from {file}: {facets.items()}")

        return facets

    @classmethod
    def decode_cmip5_netcdf(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        facets = dict()
        facets["activity"] = "CMIP5"
        facets["date"] = date
        facets["domain"] = "global"
        facets["experiment"] = data.experiment_id
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(data=data, field="frequency")
        facets["institution"] = data.institute_id
        facets["member"] = data.parent_experiment_rip
        facets["modeling_realm"] = data.modeling_realm
        facets["processing_level"] = "raw"
        facets["project"] = data.project_id
        facets["source"] = data.model_id
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderException:
            pass

        facet_schema.validate(facets)

        logging.info(f"Deciphered the following from {file}: {facets.items()}")

        return facets

    @classmethod
    def decode_cmip5_name(cls, file: Union[PathLike, str]) -> dict:
        decode_file = cls._from_filename(file=file)

        facets = dict()
        facets["activity"] = "CMIP5"
        facets["date"] = decode_file[-1]
        facets["domain"] = "global"
        facets["experiment"] = decode_file[3]
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(file=decode_file, field="frequency")
        facets["member"] = decode_file[4]
        facets["modeling_realm"] = None
        facets["processing_level"] = "raw"
        facets["source"] = decode_file[2]
        facets["variable"] = decode_file[0]
        facets["timedelta"] = cls._decode_time_info(file=decode_file, field="timedelta")
        facets["type"] = "simulation"

        try:
            logging.warning(
                "Using model to guess institute. Results may not be accurate."
            )
            facets["institution"] = CMIP5_GCM_PROVIDERS[facets["source"]]
        except KeyError:
            logging.info(f"Unable to find Institute for model: {facets['source']}")

        if "mon" in facets["frequency"]:
            facets["frequency"] = "mon"

        try:
            facets["date_start"] = date_parser(decode_file[-1])
            facets["date_end"] = date_parser(decode_file[-1], end_of_period=True)
        except DecoderException:
            pass

        facet_schema.validate(facets)

        logging.info(f"Deciphered the following from {file}: {facets.items()}")

        return facets

    @classmethod
    def decode_cordex_netcdf(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        # FIXME: What to do about our internal data that breaks all established conventions?
        facets = dict()
        facets["activity"] = "CORDEX"
        facets["date"] = date

        try:
            facets["domain"] = data.CORDEX_domain.strip()
        except AttributeError:
            try:
                facets["domain"] = data.ouranos_domain_name.strip()
            except AttributeError:
                logging.error(f"File {Path(file).name} has a nonstandard domain name.")
                raise NotImplementedError()

        facets["driving_institution"] = str(data.driving_model_id).split("-")[0]
        facets["driving_model"] = data.driving_model_id
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(data=data, field="frequency")

        if data.institute_id.strip() == "Our.":
            facets["institution"] = "Ouranos"
        else:
            facets["institution"] = data.institute_id.strip()

        facets["processing_level"] = "raw"

        if data.project_id == "":
            facets["project"] = "internal"
        else:
            facets["project"] = data.project_id

        facets["source"] = data.model_id
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except IndexError:
            pass

        try:
            facets["experiment"] = data.experiment_id.strip()
        except AttributeError:
            facets["experiment"] = data.driving_experiment_name.strip()

        try:
            facets["member"] = data.parent_experiment_rip.strip()
        except AttributeError:
            facets["member"] = data.driving_model_id.strip()

        facet_schema.validate(facets)

        logging.info(f"Deciphered the following from {file}: {facets.items()}")

        return facets

    @classmethod
    def decode_cordex_name(cls, file: Union[PathLike, str]) -> dict:
        decode_file = cls._from_filename(file=file)

        facets = dict()
        facets["activity"] = "CORDEX"
        facets["date"] = decode_file[-1]
        facets["domain"] = decode_file[1]
        facets["driving_model"] = "_".join(decode_file[2].split("-")[1:])
        facets["driving_institution"] = decode_file[2].split("-")[0]
        facets["experiment"] = decode_file[3]
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(file=decode_file, field="frequency")
        facets["institution"] = decode_file[5].split("-")[0]
        facets["member"] = decode_file[4].strip()
        facets["processing_level"] = "raw"
        facets["project"] = "CORDEX"
        facets["source"] = decode_file[5]
        facets["timedelta"] = cls._decode_time_info(file=decode_file, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = decode_file[0]

        try:
            facets["date_start"] = date_parser(decode_file[-1])
            facets["date_end"] = date_parser(decode_file[-1], end_of_period=True)
        except DecoderException:
            pass

        facet_schema.validate(facets)

        logging.info(f"Deciphered the following from {file}: {facets.items()}")

        return facets

    @classmethod
    def decode_isimip_ft_netcdf(cls, file: Union[PathLike, str]) -> dict:
        variable, date, data = cls._from_dataset(file=file)

        facets = dict()
        facets["activity"] = "ISIMP-FT"
        facets["date"] = date
        facets["co2_forcing_id"] = data.co2_forcing_id
        facets["experiment"] = data.experiment_id
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(data=data, field="frequency")
        facets["impact_model"] = data.impact_model_id
        facets["institution"] = data.institute_id
        facets["member"] = data.driving_model_ensemble_member
        facets["modeling_realm"] = data.modeling_realm
        facets["project"] = str(data.project_id)
        facets["social_forcing_id"] = data.social_forcing_id
        facets["source"] = data.model_id
        facets["timedelta"] = cls._decode_time_info(data=data, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = variable

        try:
            facets["date_start"] = date_parser(date)
            facets["date_end"] = date_parser(date, end_of_period=True)
        except DecoderException:
            pass

        facet_schema.validate(facets)

        logging.info(f"Deciphered the following from {file}: {facets.items()}")

        return facets

    @classmethod
    def decode_isimip_ft_name(cls, file: Union[PathLike, str]) -> dict:
        decode_file = cls._from_filename(file=file)

        facets = dict()
        facets["activity"] = "ISIMP-FT"
        facets["date"] = decode_file[-1]
        facets["co2_forcing_id"] = decode_file[4]
        facets["experiment"] = decode_file[2]
        facets["format"] = "netcdf"
        facets["frequency"] = cls._decode_time_info(file=decode_file, field="frequency")
        facets["impact_model_id"] = decode_file[0]
        facets["institution"] = decode_file[1].split("-")[0]
        facets["project"] = "ISIMIP-FT"
        facets["soc_forcing_id"] = decode_file[3]
        facets["source"] = "-".join(decode_file[1].split("-")[1:])
        facets["timedelta"] = cls._decode_time_info(file=decode_file, field="timedelta")
        facets["type"] = "simulation"
        facets["variable"] = decode_file[-4]

        facets["setup"] = "-".join([facets["source"], facets["experiment"]])

        if facets["co2_forcing_id"] == facets["variable"]:
            if facets["soc_forcing"] in [
                "nosoc",
                "pressoc",
                "ssp1soc",
                "ssp2",
                "ssp2soc",
                "ssp3soc",
                "ssp4soc",
                "ssp5soc",
            ]:
                facets["co2_forcing_id"] = "NAco2"
            elif facets["soc_forcing"] in ["co2", "nocco2", "pico2"]:
                facets["co2_forcing"] = facets["soc_forcing"]
                facets["soc_forcing"] = "NAsoc"

        try:
            facets["date_start"] = date_parser(decode_file[-1])
            facets["date_end"] = date_parser(decode_file[-1], end_of_period=True)
        except DecoderException:
            pass

        facet_schema.validate(facets)

        logging.info(f"Deciphered the following from {file}: {facets.items()}")

        return facets
