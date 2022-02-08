from typing import Union

import cftime
import pandas as pd
from pandas.tseries import offsets

__all__ = ["date_parser"]


def date_parser(
    date: str,
    *,
    end_of_period: bool = False,
    out_dtype: str = "str",
    strtime_format: str = "%Y-%m-%d",
) -> Union[str, pd.Timestamp]:
    """Returns a datetime from a string.

    Parameters
    ----------
    date : str
      Date to be converted
    end_of_period : bool
      If True, the date will be the end of month or year depending on what's most appropriate
    out_dtype: {"datetime", "str"}
      Returned object type.
    strtime_format: str
      If out_dtype=='str', this sets the strftime format.

    Returns
    -------
    pd.Timestamp, str
      Parsed date

    Notes
    -----
    Adapted from code written by Gabriel Rondeau-Genesse (@RondeauG)
    """

    # Formats, ordered depending on string length
    formats = {
        4: ["%Y"],
        6: ["%Y%m"],
        7: ["%Y-%m"],
        8: ["%Y%m%d"],
        10: ["%Y%m%d%H", "%Y-%m-%d"],
        12: ["%Y%m%d%H%M"],
        13: ["%Y%m-%Y%m"],
        17: ["%Y%m%d-%Y%m%d"],
        19: ["%Y-%m-%dT%H:%M:%S"],
        21: ["%Y%m%d%H-%Y%m%d%H"],
    }
    end_date_found = False

    def _parse_date(d, fmts):
        for fmt in fmts:
            try:
                s = pd.to_datetime(d, format=fmt)
                match = fmt
                break
            except Exception:  # TODO: Define this expected exception
                pass
        else:
            raise ValueError(f"Can't parse date {d} with formats {fmts}.")
        return s, match

    date_format = None
    if isinstance(date, str):
        if len(date) in [13, 17]:
            dates = date.split("-")
            if not end_of_period:
                date = dates[0]
            else:
                date = dates[1]
                end_date_found = True

        date, date_format = _parse_date(date, formats[len(date)])
    elif isinstance(date, cftime.datetime):
        for n in range(3):
            try:
                date = pd.Timestamp.fromisoformat((date - pd.Timedelta(n)).isoformat())
            except ValueError:  # We are NOT catching OutOfBoundsDatetime.
                pass
            else:
                break
        else:
            raise ValueError(
                "Unable to parse cftime date {date}, even when moving back 2 days."
            )
    elif not isinstance(date, pd.Timestamp):
        date = pd.Timestamp(date)

    if end_of_period and date_format and not end_date_found:
        if "m" not in date_format:
            date = date + pd.tseries.offsets.YearEnd(1)
        elif "d" not in date_format:
            date = date + pd.tseries.offsets.MonthEnd(1)
        # TODO: Implement sub-daily?

    if out_dtype == "str":
        return date.strftime(strtime_format)

    return date
