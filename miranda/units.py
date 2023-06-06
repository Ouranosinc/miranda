"""Special Time Units-Handling submodule."""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import xarray as xr
from xclim.core.calendar import parse_offset

KiB = int(pow(2, 10))
MiB = int(pow(2, 20))
GiB = int(pow(2, 30))


def get_time_frequency(
    d: xr.Dataset,
    expected_period: str | None = None,
    minimum_continuous_period: str = "1M",
) -> tuple[list[int | str], str]:
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
        An xarray.Dataset.
    expected_period : str
        An xarray-compatible time period (e.g. "1H", "1D", "7D", "1M", "1A").
        The time period expected of the input dataset.
        The "1M" period is specially-handled.
    minimum_continuous_period : str
        An xarray-compatible time period (e.g. "1H", "1D", "7D", "1M", "1A").
        The minimum expected granular period that data should have continuous values for.
        The "1M" period is specially-handled.

    Returns
    -------
    offset : list of int or str
        The offset a list of (multiplier, base)
    offset_meaning : str
        The offset meaning (single word)

    """
    if expected_period is not None:
        if not [expected_period.endswith(end) for end in ["H", "D", "M", "A"]]:
            raise ValueError(f"Expected period (`{expected_period}`) not supported.")

    freq = xr.infer_freq(d.time)

    # Hacky workaround for irregular Monthly data
    if freq is None or (1 < int(parse_offset(freq)[0]) < 32 and freq.endswith("D")):
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
                e_period = parse_offset(expected_period)[1]
                min_period = parse_offset(minimum_continuous_period)[1]
                collected_freqs = []
                problem_periods = []

                if e_period != "M" and min_period != "M":
                    if pd.Timedelta(expected_period) > pd.Timedelta(
                        minimum_continuous_period
                    ):
                        minimum_continuous_period = expected_period
                elif e_period == "M":
                    if pd.Timedelta(minimum_continuous_period) < pd.Timedelta(28, "D"):
                        minimum_continuous_period = expected_period

                time_periods, datasets = zip(
                    *d.time.resample(time=minimum_continuous_period)
                )

                for period, ds_part in zip(time_periods, datasets):
                    if len(ds_part) == 1:
                        logging.info(f"Skipping {str(np.datetime_as_string(period))}.")
                        # In the event that a deaccumulation/shift has created a period with one data value,
                        # we are safe in ignoring this.
                        continue

                    try:
                        f = xr.infer_freq(ds_part)
                    except ValueError as e:
                        raise ValueError(f"Issues found with {period}.") from e

                    if f is None:
                        problem_periods.append(str(np.datetime_as_string(period)))

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

    offset = [int(parse_offset(freq)[0]), parse_offset(freq)[1]]

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
