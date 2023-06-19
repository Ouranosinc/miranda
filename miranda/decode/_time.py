from __future__ import annotations

import logging
from logging import config

from pandas._libs.tslibs import NaTType  # noqa

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "DecoderError",
    "FREQUENCY_TO_POTENTIAL_TIME_UNITS",
    "TIME_UNITS_TO_FREQUENCY",
    "TIME_UNITS_TO_TIMEDELTA",
]

TIME_UNITS_TO_FREQUENCY = {
    "subhrPt": "sub-hr",
    "hourly": "1hr",
    "hours": "1hr",
    "hour": "1hr",
    "hr": "1hr",
    "1hr": "1hr",
    "3-hourly": "3hr",
    "3hr": "3hr",
    "6-hourly": "6hr",
    "6hr": "6hr",
    "Eday": "day",
    "daily": "day",
    "days": "day",
    "day": "day",
    "weekly": "sem",
    "weeks": "sem",
    "sem": "sem",
    "byweekly": "2sem",
    "2sem": "2sem",
    "monthly": "mon",
    "months": "mon",
    "month": "mon",
    "mon": "mon",
    "monC": "monC",
    "Amon": "mon",
    "Omon": "mon",
    "qtr": "3mon",
    "quarter": "3mon",
    "3mon": "3mon",
    "2qtr": "6mon",
    "semi-annual": "6mon",
    "half-yearly": "6mon",
    "6mon": "6mon",
    "yearly": "yr",
    "years": "yr",
    "annual": "yr",
    "yr": "yr",
    "yrPt": "yrPt",
    "decadal": "dec",
    "decades": "dec",
    "dec": "dec",
    "fixed": "fx",
    "fx": "fx",
}

FREQUENCY_TO_POTENTIAL_TIME_UNITS = dict()
for key, value in TIME_UNITS_TO_FREQUENCY.items():
    FREQUENCY_TO_POTENTIAL_TIME_UNITS.setdefault(value, list()).append(key)

TIME_UNITS_TO_TIMEDELTA = {
    "hourly": "1h",
    "hours": "1h",
    "hour": "1h",
    "hr": "1h",
    "1hr": "1h",
    "1hrCM": "1h",
    "1hrPt": "1h",
    "3hr": "3h",
    "3hrPt": "3h",
    "6-hourly": "6h",
    "6hr": "6h",
    "6hrPt": "6h",
    "daily": "1d",
    "day": "1d",
    "days": "1d",
    "weekly": "7d",
    "weeks": "7d",
    "sem": "7d",
    "biweekly": "14d",
    "2sem": "14d",
    "mon": "30d",
    "monC": "30d",
    "monPt": "30d",
    "Amon": "30d",
    "qtr": "90d",
    "quarter": "90d",
    "3mon": "90d",
    "2qtr": "180d",
    "6mon": "180d",
    "half-yearly": "180d",
    "semi-annual": "180d",
    "annual": "365d",
    "yearly": "365d",
    "years": "365d",
    "year": "365d",
    "yr": "365d",
    "yrPt": "365d",
}


class DecoderError(Exception):
    pass
