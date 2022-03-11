from ._institutions import *
from ._models import *
from .cmip import *
from .cmip import CMIP6_ACTIVITIES
from .cordex import *

ACTIVITIES = list()
ACTIVITIES.extend(CMIP6_ACTIVITIES)
ACTIVITIES.extend(["ERA"])


WCRP_FREQUENCIES = [
    "subhrPt",
    "1hr",
    "1hrCM",
    "1hrPt",
    "3hr",
    "3hrPt",
    "6hr",
    "6hrPt",
    "day",
    "sem",
    "mon",
    "monC",
    "monPT",
    "yr",
    "yrPt",
    "dec",
    "fx",
]
