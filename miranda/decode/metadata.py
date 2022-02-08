#!/bin/env python3
import logging
from logging import config
from os import PathLike
from pathlib import Path
from typing import List, Union

from netCDF4 import Dataset
from schema import Schema

from miranda.scripting import LOGGING_CONFIG

config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "decode_dimsvar",
    "decode_cmip5_name",
    "decode_cmip5_netcdf",
    "decode_cmip6_name",
    "decode_cmip6_netcdf",
    "decode_cordex_name",
    "decode_cordex_netcdf",
    "decode_isimip_ft_name",
    "decode_isimip_ft_netcdf",
]

schema = Schema([])


def _from_netcdf(file: Union[Path, str]) -> (str, Dataset):
    file_name = Path(file).stem
    variable_name = file_name.split("_")[0]
    data = Dataset(file)
    return variable_name, data


def _from_filename(file: Union[Path, str]) -> List[str]:
    file_name = Path(file).stem
    decode_file = file_name.split("_")
    return decode_file


CORDEX_INSTITUTES = {
    "CanRCM4": "CCCMA",
    "CRCM5": "UQAM",
    "CRCM5-UQAM": "UQAM",
    "HIRHAM5": "DMI",
    "RCA4": "SMHI",
    "RegCM4": "ISU",
    "WRF": "NCAR",
}

CMIP5_INSTITUTES = {
    "CSIRO-BOM": ["ACCESS1-0", "ACCESS1-3"],
    "BCC": ["bcc-csm1-1"],
    "CCCMA": ["CanAM4", "CanCM4", "CanESM2"],
    "CNRM-CERFACS": ["CNRM-CM5", "CNRM-CM5-2"],
    "CSIRO-QCCCE": ["CSIRO-Mk3-6-0"],
    "MOHC": ["HadCM3", "HadGEM2-ES"],
    "INM": ["inmcm4"],
    "IPSL": ["IPSL-CM5A-LR"],
    "MIROC": ["MIROC4h"],
    "MPI-M": ["MPI-ESM-LR"],
    "MRI": ["MRI-CGCM3"],
    "NCC": ["NorESM1-M"],
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


def decode_dimsvar(file: Union[Path, str]) -> dict:
    """
    see: https://gist.github.com/guziy/8543562

    Parameters
    ----------
    file: Union[Path, str]

    Returns
    -------
    dict
    """
    _, data = _from_netcdf(file=file)

    dimsvar_dict = dict()
    for var_name, varin in data.variables.items():
        dimsvar_dict[var_name] = {k: varin.getncattr(k) for k in varin.ncattrs()}

    return dimsvar_dict


def decode_era5_netcdf(file: Union[PathLike, str]) -> dict:
    pass


def decode__netcdf(file: Union[PathLike, str]) -> dict:
    pass


def decode_cmip6_netcdf(file: Union[PathLike, str]) -> dict:
    variable, data = _from_netcdf(file=file)

    facets = dict()
    facets["domain"] = "global"
    facets["experiment"] = data.experiment_id
    facets["frequency"] = data.frequency
    facets["institution"] = data.institution_id
    facets["member"] = data.variant_label
    facets["model"] = data.source_id
    facets["modeling_realm"] = data.realm
    facets["project_id"] = data.mip_era
    facets["variable"] = variable
    facets["version"] = data.version

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cmip6_name(file: Union[PathLike, str]) -> dict:
    decode_file = _from_filename(file=file)

    facets = dict()
    facets["date"] = decode_file[-1]
    facets["domain"] = "global"
    facets["experiment"] = decode_file[3]
    facets["frequency"] = decode_file[1]
    facets["grid_label"] = decode_file[5]
    facets["member"] = decode_file[4]
    facets["model"] = decode_file[2]
    facets["variable"] = decode_file[0]

    if "mon" in facets["frequency"]:
        facets["frequency"] = "mon"

    try:
        facets["institution"] = CMIP6_GCM_PROVIDERS[facets["model_id"]]
    except KeyError:
        logging.info(f"Unable to find Institute for model: {facets['model']}")

    try:
        facets["date_start"] = decode_file[6].split("-")[0]
        facets["date_end"] = decode_file[6].split("-")[1]
    except IndexError:
        pass

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cmip5_netcdf(file: Union[PathLike, str]) -> dict:
    variable, data = _from_netcdf(file=file)

    facets = dict()
    facets["domain"] = "global"
    facets["experiment"] = data.experiment_id
    facets["frequency"] = data.frequency
    facets["institution"] = data.institute_id
    facets["member"] = data.parent_experiment_rip
    facets["model"] = data.model_id
    facets["modeling_realm"] = data.modeling_realm
    facets["project"] = data.project_id
    facets["variable"] = variable

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cmip5_name(file: Union[PathLike, str]) -> dict:
    decode_file = _from_filename(file=file)

    facets = dict()
    facets["date"] = decode_file[-1]
    facets["domain"] = "global"
    facets["experiment"] = decode_file[3]
    facets["frequency"] = decode_file[1]
    facets["institution"] = CMIP5_GCM_PROVIDERS[facets["model"]]
    facets["member"] = decode_file[4]
    facets["model"] = decode_file[2]
    facets["modeling_realm"] = None
    facets["variable"] = decode_file[0]

    if "mon" in facets["frequency"]:
        facets["frequency"] = "mon"

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cordex_netcdf(file: Union[PathLike, str]) -> dict:
    variable, data = _from_netcdf(file=file)

    facets = dict()
    facets["domain"] = data.CORDEX_domain
    facets["driving_institution"] = str(data.driving_model_id).split("-")[0]
    facets["driving_model"] = data.driving_model_id
    facets["frequency"] = data.frequency
    facets["institution"] = data.institute_id
    facets["model"] = data.model_id
    facets["project"] = data.project_id
    facets["variable"] = variable

    try:
        facets["experiment"] = data.experiment_id
    except KeyError:
        facets["experiment"] = data.driving_experiment_name

    try:
        facets["member"] = data.parent_experiment_rip
    except KeyError:
        facets["member"] = data.driving_model_ensemble_member

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_cordex_name(file: Union[PathLike, str]) -> dict:
    decode_file = _from_filename(file=file)

    facets = dict()
    facets["date"] = decode_file[-1]
    facets["domain"] = decode_file[1]
    facets["driving_institution"] = "_".join(decode_file[2].split("-")[1:])
    facets["driving_institution"] = decode_file[2].split("-")[0]
    facets["experiment"] = decode_file[3]
    facets["frequency"] = decode_file[-2]
    facets["institution"] = decode_file[5].split("-")[0]
    facets["member"] = decode_file[4]
    facets["model"] = decode_file[5].split("-")[1:]
    facets["project"] = "CORDEX"
    facets["variable"] = decode_file[0]

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_isimip_ft_name(file: Union[PathLike, str]) -> dict:
    decode_file = _from_filename(file=file)

    facets = dict()
    facets["date"] = decode_file[-1]
    facets["co2_forcing_id"] = decode_file[4]
    facets["experiment"] = decode_file[2]
    facets["frequency"] = decode_file[-3]
    facets["impact_model_id"] = decode_file[0]
    facets["institution"] = decode_file[1].split("-")[0]
    facets["model"] = "-".join(decode_file[1].split("-")[1:])
    facets["project"] = "ISIMIP-FT"
    facets["setup"] = "-".join([facets["model"], facets["experiment"]])
    facets["soc_forcing_id"] = decode_file[3]
    facets["variable"] = decode_file[-4]

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

    logging.info(f"Deciphered the following from {file}: {facets.items()}")

    return facets


def decode_isimip_ft_netcdf(file: Union[PathLike, str]) -> dict:
    variable, data = _from_netcdf(file=file)

    facets = dict()
    facets["co2_forcing_id"] = data.co2_forcing_id
    facets["experiment"] = data.experiment_id
    facets["frequency"] = data.time_frequency_id
    facets["impact_model"] = data.impact_model_id
    facets["institution"] = data.institute_id
    facets["member"] = data.driving_model_ensemble_member
    facets["model"] = data.model_id
    facets["modeling_realm"] = data.modeling_realm
    facets["project"] = str(data.project_id)
    facets["social_forcing_id"] = data.social_forcing_id
    facets["variable"] = variable

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
