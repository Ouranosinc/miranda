from miranda.ecmwf._namers import ECMWF_PROJECT_NAMES  # noqa

from .cmip import CMIP5_MODELS, CMIP6_MODELS
from .cordex import CORDEX_MODELS

__all__ = ["REANALYSIS", "PROJECT_MODELS"]

REANALYSIS = list()
REANALYSIS.extend(ECMWF_PROJECT_NAMES)

PROJECT_MODELS = dict(
    CMIP5=CMIP5_MODELS, CMIP6=CMIP6_MODELS, CORDEX=CORDEX_MODELS, REANALYSIS=REANALYSIS
)

# CMIP5_INSTITUTES = {
#     "BCC": ["bcc-csm1-1"],
#     "CCCMA": ["CanAM4", "CanCM4", "CanESM2"],
#     "CNRM-CERFACS": ["CNRM-CM5", "CNRM-CM5-2"],
#     "CSIRO-BOM": ["ACCESS1-0", "ACCESS1-3"],
#     "CSIRO-QCCCE": ["CSIRO-Mk3-6-0"],
#     "ECMWF": ["ERAINT"],
#     "ICHEC": ["EC-EARTH"],
#     "INM": ["inmcm4"],
#     "IPSL": ["IPSL-CM5A-LR"],
#     "MIROC": ["MIROC4h"],
#     "MOHC": ["HadCM3", "HadGEM2-ES"],
#     "MPI-M": ["MPI-ESM-LR", "MPI-ESM-MR"],
#     "MRI": ["MRI-CGCM3"],
#     "NCC": ["NorESM1-M"],
#     "NOAA-GFDL": ["GFDL-ESM2M"],
#     "UQAM": ["GEMatm-Can", "GEMatm-MPI"],
# }
#
# CMIP6_INSTITUTES = {
#     "AWI": ["AWI-CM-1-1-MR"],
#     "BCC": ["BCC-CSM2-MR", "BCC-ESM1"],
#     "CAMS": ["CAMS-CSM1-0"],
#     "CAS": ["FGOALS-f3-L"],
#     "CCCR-IITM": ["IITM-ESM"],
#     "CCCma": ["CanESM5"],
#     "CMCC": ["CMCC-CM2-HR4", "CMCC-CM2-VHR4"],
#     "CNRM-CERFACS": ["CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1"],
#     # "DKRZ": ["MPI-ESM1-2-HR"],
#     "E3SM-Project": ["E3SM-1-0"],
#     "EC-Earth-Consortium": ["EC-Earth3", "EC-Earth3-Veg"],
#     "ECMWF": ["ECMWF-IFS-HR", "ECMWF-IFS-LR"],
#     "IPSL": ["IPSL-CM6A-ATM-HR", "IPSL-CM6A-LR"],
#     "MIROC": ["NICAM16-7S", "MIROC-ES2L", "MIROC6"],
#     "MOHC": [
#         "UKESM1-0-LL",
#         "HadGEM3-GC31-HM",
#         "HadGEM3-GC31-LL",
#         "HadGEM3-GC31-LM",
#         "HadGEM3-GC31-MM",
#     ],
#     "MPI-M": ["MPI-ESM1-2-HR", "MPI-ESM1-2-XR"],
#     "MRI": ["MRI-AGCM3-2-H", "MRI-AGCM3-2-S", "MRI-ESM2-0"],
#     "NASA-GISS": ["GISS-E2-1-G", "GISS-E2-1-H"],
#     "NCAR": ["CESM2", "CESM2-WACCM"],
#     "NOAA-GFDL": ["GFDL-AM4", "GFDL-CM4", "GFDL-CM4C192", "GFDL-ESM4", "GFDL-OM4p5B"],
#     "NUIST": ["NESM3"],
#     "SNU": ["SAM0-UNICON"],
# }
# CMIP5_GCM_PROVIDERS = {i: cat for (cat, ids) in CMIP5_INSTITUTES.items() for i in ids}
# CMIP6_GCM_PROVIDERS = {i: cat for (cat, ids) in CMIP6_INSTITUTES.items() for i in ids}
