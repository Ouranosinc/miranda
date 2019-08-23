#!/bin/env python3
import logging
from datetime import datetime as dt
from pathlib import Path

CMIP5_GCM_PROVIDERS = {
    "ACCESS1-0": "CSIRO-BOM",
    "ACCESS1-3": "CSIRO-BOM",
    "bcc-csm1-1": "BCC",
    "CanAM4": "CCCMA",
    "CanCM4": "CCCMA",
    "CanESM2": "CCCMA",
    "CNRM-CM5": "CNRM-CERFACS",
    "CNRM-CM5-2": "CNRM-CERFACS",
    "CSIRO-Mk3-6-0": "CSIRO-QCCCE",
    "HadCM3": "MOHC",
    "HadGEM2-ES": "MOHC",
    "inmcm4": "INM",
    "IPSL-CM5A-LR": "IPSL",
    "MIROC4h": "MIROC",
    "MPI-ESM-LR": "MPI-M",
    "MRI-CGCM3": "MRI",
    "NorESM1-M": "NCC",
}


def decode_cmip5(file):
    d = {}
    decode_file = Path(file).stem
    decode_file = decode_file.split("_")
    d["variable"] = decode_file[0]
    d["frequency"] = decode_file[1]
    if d["frequency"] == "Amon":
        d["frequency"] = "mon"
    d["model"] = decode_file[2]
    d["experiment"] = decode_file[3]
    d["ensemble"] = decode_file[4]
    logging.info(
        "{}: Deciphered the following from {}: {}".format(
            dt.now().strftime("%Y-%m-%d %X"), file, d.items()
        )
    )
    return d


def decode_cordex(file):
    pass


def decode_isimip_ft(file):
    pass
