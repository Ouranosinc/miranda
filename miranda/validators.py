import re

import pandas as pd
from pandas._libs.tslibs import NaTType  # noqa
from schema import Literal, Optional, Or, Regex, Schema

from miranda.metadata import CMIP6_ACTIVITIES, INSTITUTIONS, WCRP_FREQUENCIES

TYPE_NAMES = ["simulation", "reanalysis", "forecast", "gridded-obs", "station-obs"]
PROCESSING_LEVELS = ["raw", "biasadjusted"]
BASIC_DT_VALIDATION = r"\s*(?=\d{2}(?:\d{2})?)"
DATE_VALIDATION = r"^\d{4}-(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])$"

FACETS_SCHEMA = Schema(
    {
        Literal(
            "type",
            description="An Ouranos internal code used for classifying datasets.",
        ): Or(*TYPE_NAMES),
        Optional(
            Literal(
                "project",
                description="The parent project name according to the institute/authority coordinating it.",
            )
        ): str,
        Literal(
            "activity",
            description="The common climate modelling activity. "
            "Derived from 'activity_id' in WCRP-CMIP CVs",
        ): Or(*CMIP6_ACTIVITIES),
        Literal(
            "institution",
            description="The institution that created the dataset. Derived from 'institution_id' in WCRP-CMIP CVs.",
        ): str,
        "source": str,
        Optional(
            Literal(
                "driving_institution",
                description="Institute name of the global climate model driver data. "
                "Specific to regional climate model metadata.",
            )
        ): str,
        Optional(
            Literal(
                "driving_model",
                description="Model name of the global climate model driver data. "
                "Specific to regional climate model metadata.",
            )
        ): str,
        Optional(
            Literal(
                "experiment",
                description="The common experiment name. "
                "Derived from 'experiment_id' in WCRP-CMIP CVs.",
            )
        ): str,
        "frequency": Or(*WCRP_FREQUENCIES),
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
        Or("institution", "driving_institution"): Or(*INSTITUTIONS),
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
