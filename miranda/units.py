from typing import List, Optional, Tuple, Union

import pandas as pd
import xarray as xr
from xclim.core import calendar

KiB = int(pow(2, 10))
MiB = int(pow(2, 20))
GiB = int(pow(2, 30))


def get_time_frequency(
    d: xr.Dataset,
    expected_period: Optional[str] = None,
    minimum_continuous_period: str = "M",
) -> Tuple[List[Union[int, str]], str]:
    """Try to understand the Dataset frequency.

    If it can't be inferred with :py:func:`xarray.infer_freq` it tries to:
    - look for a "freq" attrs in the global or time variable attributes.
    - infer monthly frequency if all time steps are between 27 and 32 days

    In the event that an `expected_period` is supplied, special handling will be called allowing for determining data
    that may be internally discontinuous (e.g. discontinuous overall, but continuous for `minimum_continuous_period`).
    This is provided for instances where input data in a multifile dataset is sparse.

    Parameters
    ----------
    d : xr.Dataset
    expected_period : {"1H", "3H", "6H", "1D", "7D", "1M", "1Y"}
        The time period expected of the input dataset.
    minimum_continuous_period : str
        An xarray-compatible time period (e.g. "1H", "1D", "7D", "1M", "1Y")
        The minimum expected granular period that data should have continuous values for.

    Returns
    -------
    offset : List[Union[int, str]]
        The offset a list of (multiplier, base)
    offset_meaning : str
        The offset meaning (single word)
    """
    if expected_period is not None:
        if expected_period not in ["1H", "3H", "6H", "1D", "7D", "1M", "1Y"]:
            raise ValueError(f"Expected period (`{expected_period}`) not supported.")

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
            if expected_period:
                collected_freqs = []
                problem_periods = []
                format_str = {
                    "1H": "%Y-%m",
                    "3H": "%Y-%m",
                    "6H": "%Y-%m",
                    "1D": "%Y-%m",
                    "7D": "%Y-%m",
                    "1M": "%Y",
                    "1Y": "%Y",
                }[expected_period]
                if pd.Timedelta(expected_period) > pd.Timedelta(
                    minimum_continuous_period
                ):
                    minimum_continuous_period = expected_period

                time_periods, datasets = zip(
                    *d.time.resample(time=minimum_continuous_period)
                )
                for period, ds_part in zip(time_periods, datasets):
                    f = xr.infer_freq(ds_part)
                    if f is None:
                        if len(ds_part) == 1:
                            # In the event that a deaccumulation/shift has created a period with one data value,
                            # we are safe in ignoring this.
                            continue
                        else:
                            problem_periods.append(period.strftime(format_str))
                    if (
                        (d.time.diff("time") < pd.Timedelta(32, "D"))
                        & (d.time.diff("time") > pd.Timedelta(27, "D"))
                    ).all():
                        f = "1M"
                    collected_freqs.append(f)

                if problem_periods:
                    raise ValueError(
                        "Dataset contains internally discontinuous time periods: "
                        f"{' ,'.join(problem_periods)}."
                    )
                if len(set(collected_freqs)) > 1:
                    raise ValueError(
                        "Somehow, dataset contains mixed frequencies: "
                        f"{' ,'.join(collected_freqs)}."
                    )
                freq = set(collected_freqs).pop()
            else:
                raise ValueError("Dataset time component may be discontinuous.")

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
