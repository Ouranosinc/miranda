import re

import pandas as pd
from pandas._libs.tslibs import NaTType  # noqa
from schema import Literal, Optional, Or, Regex, Schema

TYPE_NAMES = ["simulation", "reanalysis", "forecast", "gridded-obs", "station-obs"]
PROCESSING_LEVELS = ["raw", "biasadjusted"]
BASIC_DT_VALIDATION = r"\s*(?=\d{2}(?:\d{2})?)"
DATE_VALIDATION = r"^\d{4}-(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])$"
CMIP6_FREQUENCIES = [
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


FACETS_SCHEMA = Schema(
    {
        Literal(
            "type",
            description="An Ouranos internal code used for classifying datasets.",
        ): Or(TYPE_NAMES),
        Optional(
            Literal(
                "project",
                description="The parent project name according to the institute/authority coordinating it.",
            )
        ): str,
        "activity": str,
        Literal(
            "institution",
            description="The institution that created the dataset. Typically derived from 'institution_id'.",
        ): str,
        "source": str,
        Optional("driving_institution"): str,
        Optional("driving_model"): str,
        Optional("experiment"): str,
        "frequency": Or(*CMIP6_FREQUENCIES),
        "domain": str,
        Optional("member"): str,
        Optional("variable"): str,
        Optional("timedelta"): Or(pd.Timedelta, NaTType, "NaT"),
        Optional("date"): Or(Regex(BASIC_DT_VALIDATION, flags=re.I), "fx"),
        Optional("date_start"): Or(Regex(DATE_VALIDATION, flags=re.I), NaTType, "NaT"),
        Optional("date_end"): Or(Regex(DATE_VALIDATION, flags=re.I), NaTType, "NaT"),
        Optional("processing_level"): Or(*PROCESSING_LEVELS),
        "format": Or("netcdf", "zarr"),
        Optional("version"): str,
    },
    ignore_extra_keys=True,
)


STATION_OBS_SCHEMA = Schema(
    {
        "type": "station-obs",
        "project": str,
        "institution": str,
        "source": str,
        "frequency": str,
        "variable": str,
        "version": str,
    },
    ignore_extra_keys=True,
)


GRIDDED_SCHEMA = Schema(
    {
        "type": Or("forecast", "gridded-obs", "reanalysis"),
        "project": str,
        "institution": str,
        "source": str,
        "domain": str,
        "frequency": str,
        "variable": str,
    },
    ignore_extra_keys=True,
)

SIMULATION_SCHEMA = Schema(
    {
        "type": "simulation",
        "processing_level": Or(*PROCESSING_LEVELS),
        "project": str,
        "institution": str,
        "source": str,
        "domain": str,
        Or("driving_model", "member"): str,
        "experiment": str,
        "frequency": str,
        "member": str,
        "variable": str,
    },
    ignore_extra_keys=True,
)
