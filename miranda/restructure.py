#!/bin/env python3
import logging
from collections import defaultdict
from datetime import datetime as dt
from pathlib import Path

from netCDF4 import Dataset


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


def decode_cmip5_netcdf(file: str or Path) -> dict:
    facets = defaultdict(str)
    decode_file = Path(file).stem
    decode_file = str(decode_file).split("_")
    data = Dataset(file)

    facets["variable"] = decode_file[0]
    facets["project"] = data.project_id
    facets["institute"] = data.institute_id
    facets["frequency"] = data.frequency
    facets["model"] = data.model_id
    facets["experiment"] = data.experiment_id
    facets["ensemble"] = data.parent_experiment_rip
    facets["realm"] = data.modeling_realm

    logging.info(
        "{}: Deciphered the following from {}: {}".format(
            dt.now().strftime("%Y-%m-%d %X"), file, facets.items()
        )
    )

    return facets


def decode_cmip5_name(file: str or Path) -> dict:
    facets = dict()
    decode_file = Path(file).stem
    decode_file = str(decode_file).split("_")

    facets["variable"] = decode_file[0]
    facets["frequency"] = decode_file[1]
    if facets["frequency"] == "Amon":
        facets["frequency"] = "mon"
    facets["model"] = decode_file[2]
    facets["institute"] = CMIP5_GCM_PROVIDERS[facets["model"]]
    facets["realm"] = None
    facets["experiment"] = decode_file[3]
    facets["ensemble"] = decode_file[4]

    logging.info(
        "{}: Deciphered the following from {}: {}".format(
            dt.now().strftime("%Y-%m-%d %X"), file, facets.items()
        )
    )

    return facets


def decode_cmip6_netcdf(file: str or Path) -> dict:
    pass


def decode_cmip6_name(file: str or Path) -> dict:
    pass


def decode_cordex_netcdf(file: str or Path) -> dict:
    facets = defaultdict(str)
    decode_file = Path(file).stem
    decode_file = str(decode_file).split("_")
    data = Dataset(file)

    facets["variable"] = decode_file[0]
    facets["project"] = data.project_id
    facets["institute"] = data.institute_id
    facets["model"] = data.model_id
    facets["CORDEX_domain"] = data.CORDEX_domain
    facets["frequency"] = data.frequency
    facets["driver"] = data.driving_model_id
    facets["experiment"] = data.experiment_id
    facets["ensemble"] = data.parent_experiment_rip
    facets["realm"] = data.modeling_realm

    logging.info(
        "{}: Deciphered the following from {}: {}".format(
            dt.now().strftime("%Y-%m-%d %X"), file, facets.items()
        )
    )

    return facets


def decode_cordex_name(file: str or Path) -> dict:
    facets = dict()
    decode_file = Path(file).stem
    decode_file = str(decode_file).split("_")

    facets["institute"] = decode_file[5].split("-")[0]
    facets["model"] = decode_file[5].split("-")[1:]
    facets["experiment"] = "_".join(decode_file[1:4])
    facets["time_frequency"] = decode_file[-2]
    facets["realm"] = None
    facets["ensemble"] = decode_file[4]
    facets["variable"] = decode_file[0]

    logging.info(
        "{}: Deciphered the following from {}: {}".format(
            dt.now().strftime("%Y-%m-%d %X"), file, facets.items()
        )
    )

    return facets


def decode_isimip_ft_name(file: str or Path) -> dict:
    facets = dict()
    f = Path(file)
    segm = str(f.stem).split("_")

    facets["variable"] = segm[-4]
    facets["project"] = "ISIMIP-FT"
    facets["institute"] = segm[1].split("-")[0]
    facets["model"] = "-".join(segm[1].split("-")[1:])
    facets["experiment"] = segm[2]

    facets["setup"] = "-".join([facets["model"], facets["experiment"]])

    facets["time_frequency"] = segm[-3]
    facets["impact_model"] = segm[0]
    facets["soc_forcing"] = segm[3]
    facets["co2_forcing"] = segm[4]

    if facets["co2_forcing"] == facets["variable"]:
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
            facets["co2_forcing"] = "NAco2"
        elif facets["soc_forcing"] in ["co2", "nocco2", "pico2"]:
            facets["co2_forcing"] = facets["soc_forcing"]
            facets["soc_forcing"] = "NAsoc"

    logging.info(
        "{}: Deciphered the following from {}: {}".format(
            dt.now().strftime("%Y-%m-%d %X"), file, facets.items()
        )
    )

    return facets


def decode_isimip_ft_netcdf(file):
    facets = defaultdict(str)
    decode_file = Path(file).stem
    decode_file = str(decode_file).split("_")
    data = Dataset(file)

    facets["variable"] = decode_file[0]
    facets["project"] = str(data.project_id).upper()
    facets["institute"] = data.institute_id
    facets["frequency"] = data.time_frequency_id
    facets["model"] = data.model_id
    facets["impact_model"] = data.impact_model_id
    facets["social_forcing"] = data.social_forcing_id
    facets["co2_forcing"] = data.co2_forcing_id
    facets["experiment"] = data.experiment_id
    facets["ensemble"] = data.driving_model_ensemble_member
    facets["realm"] = data.modeling_realm

    logging.info(
        "{}: Deciphered the following from {}: {}".format(
            dt.now().strftime("%Y-%m-%d %X"), file, facets.items()
        )
    )

    return facets
