from typing import List, Tuple, Union

import pandas as pd
import xarray as xr
from xclim.core import calendar

KiB = int(pow(2, 10))
MiB = int(pow(2, 20))
GiB = int(pow(2, 30))


def get_time_frequency(d: xr.Dataset) -> Tuple[List[Union[int, str]], str]:
    """Try to understand the Dataset frequency.

    If it can't be inferred with :py:func:`xarray.infer_freq` it tries to:
    - look for a "freq" attrs in the global or time variable attributes.
    - infer monthly frequency if all time steps are between 27 and 32 days

    returns the offset a list of (multiplier, base) and its meaning (single word)
    """
    freq = xr.infer_freq(d.time)

    # Hacky workaround for irregular Monthly data
    if freq is None or (
        1 < int(calendar.parse_offset(freq)[0]) < 32 and freq.endswith("D")
    ):
        if "freq" in d.attrs:
            freq = d.attrs["freq"]
        elif "freq" in d.time.attrs:
            freq = d.time.attrs["freq"]
        elif (
            (d.time.diff("time") < pd.Timedelta(32, "D"))
            & (d.time.diff("time") > pd.Timedelta(27, "D"))
        ).all():
            freq = "1M"
        else:
            raise TypeError("Dataset time component may be discontinuous.")

    offset = [int(calendar.parse_offset(freq)[0]), calendar.parse_offset(freq)[1]]

    time_units = {
        "s": "second",
        "T": "minute",
        "h": "hour",
        "D": "day",
        "M": "month",
        "W": "week",
        "A": "year",
    }
    if offset[1] in ["S", "H"]:
        offset[1] = offset[1].lower()
    offset_meaning = time_units[offset[1]]
    return offset, offset_meaning
