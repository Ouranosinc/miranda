#!/bin/env python3
import logging
import re
from logging import config
from os import PathLike
from pathlib import Path
from typing import List, Union

import netCDF4 as nc
import pandas as pd
import schema
import zarr
from pandas._libs.tslibs import NaTType  # noqa

from miranda.scripting import LOGGING_CONFIG

from ._utils import date_parser

config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "CMIP5_GCM_PROVIDERS",
    "CMIP5_INSTITUTES",
    "CMIP6_GCM_PROVIDERS",
    "CMIP6_INSTITUTES",
    "decode_eccc_obs",
    "decode_ahccd_obs",
    "decode_melcc_obs",
    "decode_generic_reanalysis",
    "decode_era5",
    "decode_cmip5_name",
    "decode_cmip5_netcdf",
    "decode_cmip6_name",
    "decode_cmip6_netcdf",
    "decode_cordex_name",
    "decode_cordex_netcdf",
    "decode_primary_variable",
    "decode_isimip_ft_name",
    "decode_isimip_ft_netcdf",
]

BASIC_DT_VALIDATION = r"\s*(?=\d{2}(?:\d{2})?)"
DATE_VALIDATION = r"^\d{4}-(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])$"
FREQUENCY_TO_TIMEDELTA = {
    "hourly": "1h",
    "6-hourly": "6h",
    "daily": "1d",
    "day": "1d",
    "weekly": "7d",
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
            str, lambda f: f in ["1hr", "3hr", "6hr", "day", "dec", "mon", "yr", "fx"]
        ),
        "domain": str,
        schema.Optional("member"): str,
        schema.Optional("variable"): str,
        schema.Optional("timedelta"): schema.Or(pd.Timedelta, NaTType),
        schema.Optional("date"): schema.Or(
            schema.Regex(BASIC_DT_VALIDATION, flags=re.I), NaTType
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
    },
    ignore_extra_keys=True,
)

CORDEX_MODELS = [
    "ALADIN52",
    "ALADIN53",
    "ALADIN63",
    "ALADIN64",
    "ALARO-0",
    "BOM-SDM",
    "CCAM",
    "CCAM-1704",
    "CCAM-2008",
    "CCLM-0-9",
    "CCLM4-21-2",
    "CCLM4-8-17",
    "CCLM4-8-17-CLM3-5",
    "CCLM5-0-15",
    "CCLM5-0-2",
    "CCLM5-0-6",
    "CCLM5-0-9",
    "COSMO-crCLIM-v1-1",
    "CRCM5",
    "CRCM5-SN",
    "CanRCM4",
    "Eta",
    "HIRHAM5",
    "HadGEM3-RA",
    "HadREM3-GA7-05",
    "HadRM3P",
    "MAR311",
    "MAR36",
    "RA",
    "RACMO21P",
    "RACMO22E",
    "RACMO22T",
    "RCA4",
    "RCA4-SN",
    "REMO2009",
    "REMO2015",
    "RRCM",
    "RegCM4",
    "RegCM4-0",
    "RegCM4-2",
    "RegCM4-3",
    "RegCM4-4",
    "RegCM4-6",
    "RegCM4-7",
    "SNURCM",
    "VRF370",
    "WRF",
    "WRF331",
    "WRF331F",
    "WRF331G",
    "WRF341E",
    "WRF341I",
    "WRF351",
    "WRF360J",
    "WRF360K",
    "WRF360L",
    "WRF361H",
    "WRF381P",
]

CMIP5_INSTITUTES = {
    "BCC": ["bcc-csm1-1"],
    "CCCMA": ["CanAM4", "CanCM4", "CanESM2"],
    "CNRM-CERFACS": ["CNRM-CM5", "CNRM-CM5-2"],
    "CSIRO-BOM": ["ACCESS1-0", "ACCESS1-3"],
    "CSIRO-QCCCE": ["CSIRO-Mk3-6-0"],
    "ECMWF": ["ERAINT"],
    "ICHEC": ["EC-EARTH"],
    "INM": ["inmcm4"],
    "IPSL": ["IPSL-CM5A-LR"],
    "MIROC": ["MIROC4h"],
    "MOHC": ["HadCM3", "HadGEM2-ES"],
    "MPI-M": ["MPI-ESM-LR", "MPI-ESM-MR"],
    "MRI": ["MRI-CGCM3"],
    "NCC": ["NorESM1-M"],
    "NOAA-GFDL": ["GFDL-ESM2M"],
    "UQAM": ["GEMatm-Can", "GEMatm-MPI"],
}

CMIP6_INSTITUTES = {
    "AWI": ["AWI-CM-1-1-MR"],
    "BCC": ["BCC-CSM2-MR", "BCC-ESM1"],
    "CAMS": ["CAMS-CSM1-0"],
    "CAS": ["FGOALS-f3-L"],
    "CCCR-IITM": ["IITM-ESM"],
    "CCCma": ["CanESM5"],
    "CMCC": ["CMCC-CM2-HR4", "CMCC-CM2-VHR4"],
    "CNRM-CERFACS": ["CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1"],
    # "DKRZ": ["MPI-ESM1-2-HR"],
    "E3SM-Project": ["E3SM-1-0"],
    "EC-Earth-Consortium": ["EC-Earth3", "EC-Earth3-Veg"],
    "ECMWF": ["ECMWF-IFS-HR", "ECMWF-IFS-LR"],
    "IPSL": ["IPSL-CM6A-ATM-HR", "IPSL-CM6A-LR"],
    "MIROC": ["NICAM16-7S", "MIROC-ES2L", "MIROC6"],
    "MOHC": [
        "UKESM1-0-LL",
        "HadGEM3-GC31-HM",
        "HadGEM3-GC31-LL",
        "HadGEM3-GC31-LM",
        "HadGEM3-GC31-MM",
    ],
    "MPI-M": ["MPI-ESM1-2-HR", "MPI-ESM1-2-XR"],
    "MRI": ["MRI-AGCM3-2-H", "MRI-AGCM3-2-S", "MRI-ESM2-0"],
    "NASA-GISS": ["GISS-E2-1-G", "GISS-E2-1-H"],
    "NCAR": ["CESM2", "CESM2-WACCM"],
    "NOAA-GFDL": ["GFDL-AM4", "GFDL-CM4", "GFDL-CM4C192", "GFDL-ESM4", "GFDL-OM4p5B"],
    "NUIST": ["NESM3"],
    "SNU": ["SAM0-UNICON"],
}

CMIP5_GCM_PROVIDERS = {i: cat for (cat, ids) in CMIP5_INSTITUTES.items() for i in ids}
CMIP6_GCM_PROVIDERS = {i: cat for (cat, ids) in CMIP6_INSTITUTES.items() for i in ids}


def _from_dataset(file: Union[Path, str]) -> (str, str, nc.Dataset):
    file_name = Path(file).stem

    variable_name = decode_primary_variable(file)
    variable_date = file_name.split("_")[-1]
    data = nc.Dataset(file)
    return variable_name, variable_date, data


def _from_filename(file: Union[Path, str]) -> List[str]:
    file_name = Path(file).stem
    decode_file = file_name.split("_")
    return decode_file


def decode_primary_variable(file: Union[Path, str]) -> str:
    """Attempts to find the primary variable of a netCDF

    Parameters
    ----------
    file: Union[Path, str]

    Returns
    -------
    str
    """

    dimsvar_dict = dict()
    coords = ("time", "lat", "lon")

    if file.is_file() and file.suffix in ["nc", "nc4"]:
        data = nc.Dataset(file, mode="r")
        for var_name, var_attrs in data.variables.items():
            dimsvar_dict[var_name] = {
                k: var_attrs.getncattr(k) for k in var_attrs.ncattrs()
            }
        for k in dimsvar_dict.keys():
            if not str(k).startswith(coords):
                return str(k)

    elif file.is_dir() and file.suffix == "zarr":
        data = zarr.open(file, mode="r")
        for k in data.array_keys():
            if not str(k).startswith(coords):
                return str(k)


def decode_era5(file: Union[PathLike, str]) -> dict:
    pass


def decode_generic_reanalysis(file: Union[PathLike, str]) -> dict:
    pass


def decode_eccc_obs(file: Union[PathLike, str]) -> dict:
    pass


def decode_ahccd_obs(file: Union[PathLike, str]) -> dict:
    pass


def decode_melcc_obs(file: Union[PathLike, str]) -> dict:
    pass


def decode_cmip6_netcdf(file: Union[PathLike, str]) -> dict:
    variable, date, data = _from_dataset(file=file)

    facets = dict()
    facets["activity"] = data.activity_id
    facets["date"] = date
    facets["domain"] = "global"
    facets["experiment"] = data.experiment_id
    facets["format"] = "netcdf"
    facets["frequency"] = data.frequency
    facets["institution"] = data.institution_id
    facets["member"] = data.variant_label
    facets["modeling_realm"] = data.realm
    facets["processing_level"] = "raw"
    facets["project"] = data.project
    facets["source"] = data.source_id
    facets["timedelta"] = pd.to_timedelta(FREQUENCY_TO_TIMEDELTA[data.frequency])
    facets["type"] = "simulation"
    facets["variable"] = variable
    facets["version"] = data.version

    try:
        facets["date_start"] = date_parser(date)
        facets["date_end"] = date_parser(date, end_of_period=True)
    except IndexError:
        pass

    facet_schema.validate(facets)

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cmip6_name(file: Union[PathLike, str]) -> dict:
    decode_file = _from_filename(file=file)

    facets = dict()
    facets["activity"] = "CMIP6"
    facets["date"] = decode_file[-1]
    facets["domain"] = "global"
    facets["experiment"] = decode_file[3]
    facets["format"] = "netcdf"
    facets["frequency"] = decode_file[1]
    facets["grid_label"] = decode_file[5]
    facets["member"] = decode_file[4]
    facets["processing_level"] = "raw"
    facets["project"] = "CMIP6"
    facets["source"] = decode_file[2]
    facets["timedelta"] = pd.to_timedelta(FREQUENCY_TO_TIMEDELTA[decode_file[1]])
    facets["type"] = "simulation"
    facets["variable"] = decode_file[0]

    if "mon" in facets["frequency"]:
        facets["frequency"] = "mon"

    try:
        facets["institution"] = CMIP6_GCM_PROVIDERS[facets["source"]]
    except KeyError:
        logging.info(f"Unable to find Institute for model: {facets['source']}")

    try:
        facets["date_start"] = date_parser(decode_file[-1])
        facets["date_end"] = date_parser(decode_file[-1], end_of_period=True)
    except IndexError:
        pass

    facet_schema.validate(facets)

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cmip5_netcdf(file: Union[PathLike, str]) -> dict:
    variable, date, data = _from_dataset(file=file)

    facets = dict()
    facets["activity"] = "CMIP5"
    facets["date"] = date
    facets["domain"] = "global"
    facets["experiment"] = data.experiment_id
    facets["format"] = "netcdf"
    facets["frequency"] = data.frequency
    facets["institution"] = data.institute_id
    facets["member"] = data.parent_experiment_rip
    facets["modeling_realm"] = data.modeling_realm
    facets["processing_level"] = "raw"
    facets["project"] = data.project_id
    facets["source"] = data.model_id
    facets["timedelta"] = pd.to_timedelta(FREQUENCY_TO_TIMEDELTA[data.frequency])
    facets["type"] = "simulation"
    facets["variable"] = variable

    try:
        facets["date_start"] = date_parser(date)
        facets["date_end"] = date_parser(date, end_of_period=True)
    except IndexError:
        pass

    facet_schema.validate(facets)

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cmip5_name(file: Union[PathLike, str]) -> dict:
    decode_file = _from_filename(file=file)

    facets = dict()
    facets["activity"] = "CMIP5"
    facets["date"] = decode_file[-1]
    facets["domain"] = "global"
    facets["experiment"] = decode_file[3]
    facets["format"] = "netcdf"
    facets["frequency"] = decode_file[-2]
    facets["member"] = decode_file[4]
    facets["modeling_realm"] = None
    facets["processing_level"] = "raw"
    facets["source"] = decode_file[2]
    facets["variable"] = decode_file[0]
    facets["timedelta"] = pd.to_timedelta(FREQUENCY_TO_TIMEDELTA[decode_file[-2]])
    facets["type"] = "simulation"

    facets["institution"] = CMIP5_GCM_PROVIDERS[facets["source"]]

    if "mon" in facets["frequency"]:
        facets["frequency"] = "mon"

    try:
        facets["date_start"] = date_parser(decode_file[-1])
        facets["date_end"] = date_parser(decode_file[-1], end_of_period=True)
    except IndexError:
        pass

    facet_schema.validate(facets)

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cordex_netcdf(file: Union[PathLike, str]) -> dict:
    variable, date, data = _from_dataset(file=file)

    facets = dict()
    facets["activity"] = "CORDEX"
    facets["date"] = date
    facets["domain"] = data.CORDEX_domain
    facets["driving_institution"] = str(data.driving_model_id).split("-")[0]
    facets["driving_model"] = data.driving_model_id
    facets["format"] = "netcdf"
    facets["frequency"] = data.frequency
    facets["institution"] = data.institute_id
    facets["processing_level"] = "raw"
    facets["project"] = data.project_id
    facets["source"] = data.model_id
    facets["timedelta"] = pd.to_timedelta(FREQUENCY_TO_TIMEDELTA[data.frequency])
    facets["type"] = "simulation"
    facets["variable"] = variable

    try:
        facets["date_start"] = date_parser(date)
        facets["date_end"] = date_parser(date, end_of_period=True)
    except IndexError:
        pass

    try:
        facets["experiment"] = data.experiment_id
    except KeyError:
        facets["experiment"] = data.driving_experiment_name

    try:
        facets["member"] = data.parent_experiment_rip
    except KeyError:
        facets["member"] = data.driving_model_ensemble_member

    facet_schema.validate(facets)

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cordex_name(file: Union[PathLike, str]) -> dict:
    decode_file = _from_filename(file=file)

    facets = dict()
    facets["activity"] = "CORDEX"
    facets["date"] = decode_file[-1]
    facets["domain"] = decode_file[1]
    facets["driving_model"] = "_".join(decode_file[2].split("-")[1:])
    facets["driving_institution"] = decode_file[2].split("-")[0]
    facets["experiment"] = decode_file[3]
    facets["format"] = "netcdf"
    facets["frequency"] = decode_file[-2]
    facets["institution"] = decode_file[5].split("-")[0]
    facets["member"] = decode_file[4]
    facets["processing_level"] = "raw"
    facets["project"] = "CORDEX"
    facets["source"] = decode_file[5]
    facets["timedelta"] = pd.to_timedelta(FREQUENCY_TO_TIMEDELTA[decode_file[-2]])
    facets["type"] = "simulation"
    facets["variable"] = decode_file[0]

    try:
        facets["date_start"] = date_parser(decode_file[-1])
        facets["date_end"] = date_parser(decode_file[-1], end_of_period=True)
    except IndexError:
        pass

    facet_schema.validate(facets)

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_isimip_ft_netcdf(file: Union[PathLike, str]) -> dict:
    variable, date, data = _from_dataset(file=file)

    facets = dict()
    facets["activity"] = "ISIMP-FT"
    facets["date"] = date
    facets["co2_forcing_id"] = data.co2_forcing_id
    facets["experiment"] = data.experiment_id
    facets["format"] = "netcdf"
    facets["frequency"] = data.time_frequency_id
    facets["impact_model"] = data.impact_model_id
    facets["institution"] = data.institute_id
    facets["member"] = data.driving_model_ensemble_member
    facets["modeling_realm"] = data.modeling_realm
    facets["project"] = str(data.project_id)
    facets["social_forcing_id"] = data.social_forcing_id
    facets["source"] = data.model_id
    facets["timedelta"] = pd.to_timedelta(FREQUENCY_TO_TIMEDELTA[data.frequency])
    facets["type"] = "simulation"
    facets["variable"] = variable

    try:
        facets["date_start"] = date_parser(date)
        facets["date_end"] = date_parser(date, end_of_period=True)
    except IndexError:
        pass

    facet_schema.validate(facets)

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_isimip_ft_name(file: Union[PathLike, str]) -> dict:
    decode_file = _from_filename(file=file)

    facets = dict()
    facets["activity"] = "ISIMP-FT"
    facets["date"] = decode_file[-1]
    facets["co2_forcing_id"] = decode_file[4]
    facets["experiment"] = decode_file[2]
    facets["format"] = "netcdf"
    facets["frequency"] = decode_file[-3]
    facets["impact_model_id"] = decode_file[0]
    facets["institution"] = decode_file[1].split("-")[0]
    facets["project"] = "ISIMIP-FT"
    facets["soc_forcing_id"] = decode_file[3]
    facets["source"] = "-".join(decode_file[1].split("-")[1:])
    facets["timedelta"] = pd.to_timedelta(FREQUENCY_TO_TIMEDELTA[decode_file[-3]])
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
    except IndexError:
        pass

    facet_schema.validate(facets)

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


# Decoder class idea
class DecoderException(Exception):
    pass


class Decoder:
    pass


class CORDEXDecoder(Decoder):
    pass


class CMIP5Decoder(Decoder):
    pass


class CMIP6Decoder(Decoder):
    pass


class ISIMIPDecorder(Decoder):
    pass
