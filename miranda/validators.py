"""Data Validation module."""
from __future__ import annotations

import re
import typing

import pandas as pd
from pandas._libs.tslibs import NaTType  # noqa
from schema import Literal, Optional, Or, Regex, Schema

from .cv import VALIDATION_ENABLED

__all__ = ["url_validate"]

if VALIDATION_ENABLED:
    from .cv import (
        ACTIVITIES,
        BIAS_ADJUST_INSTITUTIONS,
        DRIVING_MODELS,
        INSTITUTIONS,
        WCRP_FREQUENCIES,
    )

    __all__ = ["validation_schemas", "url_validate"]

    TYPE_NAMES = [
        "simulation",
        "reconstruction",
        "forecast",
        "gridded-obs",
        "station-obs",
    ]
    PROCESSING_LEVELS = ["raw", "biasadjusted", "extracted", "regridded"]
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
                    "activity",
                    description="The common climate modelling activity. "
                    "Derived from 'activity_id' in WCRP-CMIP CVs",
                )
            ): Or(*ACTIVITIES),
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
            Optional("date_start"): Or(
                Regex(DATE_VALIDATION, flags=re.I), NaTType, "NaT"
            ),
            Optional("date_end"): Or(
                Regex(DATE_VALIDATION, flags=re.I), NaTType, "NaT"
            ),
            Optional("processing_level"): Or(*PROCESSING_LEVELS),
            "format": Or("netcdf", "zarr"),
            Optional("version"): str,
        },
        ignore_extra_keys=True,
    )

    STATION_OBS_SCHEMA = Schema(
        {
            "type": "station-obs",
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
            "type": Or("forecast", "gridded-obs", "reconstruction"),
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
            "activity": str,
            "mip_era": str,
            "institution": Or(*INSTITUTIONS),
            Optional("driving_institution"): Or(*INSTITUTIONS),
            "source": str,
            "domain": str,
            Optional("driving_model"): Or(*DRIVING_MODELS),
            Optional("bias_adjust_institution"): Or(*BIAS_ADJUST_INSTITUTIONS),
            Optional("bias_adjust_project"): str,
            "experiment": str,
            "member": str,
            "frequency": str,
            "variable": str,
        },
        ignore_extra_keys=True,
    )

    validation_schemas = dict()
    validation_schemas["simulation"] = SIMULATION_SCHEMA
    validation_schemas["station-obs"] = STATION_OBS_SCHEMA
    validation_schemas.update(
        {
            data_type: GRIDDED_SCHEMA
            for data_type in ["forecast", "gridded-obs", "reconstruction"]
        }
    )


def url_validate(target: str) -> typing.Match[str] | None:
    """Validate whether a supplied URL is reliably written.

    Parameters
    ----------
    target : str

    References
    ----------
    https://stackoverflow.com/a/7160778/7322852
    """
    url_regex = re.compile(
        r"^(?:http|ftp)s?://"  # http:// or https://
        # domain...
        r"(?:(?:[A-Z\d](?:[A-Z\d-]{0,61}[A-Z\d])?\.)+(?:[A-Z]{2,6}\.?|[A-Z\d-]{2,}\.?)|"
        r"localhost|"  # localhost...
        r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"  # ...or ip
        r"(?::\d+)?"  # optional port
        r"(?:/?|[/?]\S+)$",
        re.IGNORECASE,
    )
    return re.match(url_regex, target)
