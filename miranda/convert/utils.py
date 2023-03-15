import hashlib
import logging.config
import os
import re
from pathlib import Path
from typing import Dict, Optional, Union

import cftime
import pandas as pd
import xarray as xr
from pandas._libs import NaTType
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "daily_aggregation",
    "find_version_hash",
]


def daily_aggregation(
    ds: xr.Dataset, keys_only: bool = False
) -> Dict[str, Optional[xr.Dataset]]:
    logging.info("Creating daily upscaled climate variables.")

    daily_dataset = dict()
    for variable in ds.data_vars:
        if variable in ["tas", "tdps"]:
            # Some looping to deal with memory consumption issues
            for v, func in {
                f"{variable}max": "max",
                f"{variable}min": "min",
                f"{variable}": "mean",
            }.items():
                if not keys_only:
                    ds_out = xr.Dataset()
                    ds_out.attrs = ds.attrs.copy()
                    ds_out.attrs["frequency"] = "day"

                    method = f"time: {func}{'imum' if func != 'mean' else ''} (interval: 1 day)"
                    ds_out.attrs["cell_methods"] = method

                    if v == "tas" and not hasattr(ds, "tas"):
                        ds_out[v] = tas(tasmax=ds.tasmax, tasmin=ds.tasmin)
                    else:
                        # Thanks for the help, xclim contributors
                        r = ds[variable].resample(time="D")
                        ds_out[v] = getattr(r, func)(dim="time", keep_attrs=True)

                    daily_dataset[v] = ds_out
                    del ds_out
                else:
                    daily_dataset[v] = None

        elif variable in [
            "evspsblpot",
            "hfls",
            "hfss",
            "hur",
            "hus",
            "pr",
            "prsn",
            "ps",
            "psl",
            "rsds",
            "rss",
            "rlds",
            "rls",
            "snd",
            "snr",
            "snw",
            "swe",
        ]:
            if not keys_only:
                ds_out = xr.Dataset()
                ds_out.attrs = ds.attrs.copy()
                ds_out.attrs["frequency"] = "day"
                ds_out.attrs["cell_methods"] = "time: mean (interval: 1 day)"
                logging.info(f"Converting {variable} to daily time step (daily mean).")
                ds_out[variable] = (
                    ds[variable].resample(time="D").mean(dim="time", keep_attrs=True)
                )

                daily_dataset[variable] = ds_out
                del ds_out
            else:
                daily_dataset[variable] = None
        else:
            continue

    return daily_dataset


def find_version_hash(file: Union[os.PathLike, str]) -> Dict:
    """

    Parameters
    ----------
    file : Path or str

    Returns
    -------
    dict
    """

    def _get_hash(f):
        hash_sha256_writer = hashlib.sha256()
        with open(f, "rb") as f_opened:
            hash_sha256_writer.update(f_opened.read())
        sha256sum = hash_sha256_writer.hexdigest()
        logging.info(f"Calculated sha256sum (starting: {sha256sum[:6]})")
        del hash_sha256_writer
        return sha256sum

    version_info = dict()
    possible_version = Path(file).parent.name
    if re.match(r"^v\d+", possible_version, re.IGNORECASE):
        version_info["version"] = Path(file).parent.name
        version_info["sha256sum"] = _get_hash(file)

    else:
        file_identity = str(Path(file).name).split(".")[0]
        possible_version_signature = Path(file).parent.glob(f"{file_identity}.*")
        for sig in possible_version_signature:
            found_version = re.search(r"\.(v\d+.+)$", sig.name, re.IGNORECASE)
            if found_version:
                try:
                    version_info["version"] = found_version.group()
                    version_info["sha256sum"] = int(sig.open().read())
                except ValueError:
                    continue
                break
        else:
            version_info["version"] = "vNotFound"
            version_info["sha256sum"] = _get_hash(file)

    return version_info


def date_parser(
    date: str,
    *,
    end_of_period: bool = False,
    output_type: str = "str",
    strftime_format: str = "%Y-%m-%d",
) -> Union[str, pd.Timestamp, NaTType]:
    """Parses datetime objects from a string representation of a date or both a start and end date.

    Parameters
    ----------
    date : str
        Date to be converted.
    end_of_period : bool
        If True, the date will be the end of month or year depending on what's most appropriate.
    output_type : {"datetime", "str"}
        Desired returned object type.
    strftime_format : str
        If output_type=='str', this sets the strftime format.

    Returns
    -------
    pd.Timestamp or str or pd.NaT
        Parsed date.

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
        25: ["%Y%m%d%H%M-%Y%m%d%H%M"],
    }
    end_date_found = False

    def _parse_date(d, fmts):
        for fmt in fmts:
            try:
                s = pd.to_datetime(d, format=fmt)
                match = fmt
                break
            except ValueError:
                pass
        else:
            raise ValueError(
                f"Can't parse date {d} with supported formats: [{', '.join(fmts)}]."
            )
        return s, match

    date_format = None
    if isinstance(date, str):
        if len(date) in [13, 17, 21, 25]:
            dates = date.split("-")
            if not end_of_period:
                date = dates[0]
            else:
                date = dates[1]
                end_date_found = True

        try:
            possible_formats = formats[len(date)]
            date, date_format = _parse_date(date, possible_formats)
        except KeyError:
            # Return NaT for fixed/missing/ill-formatted date strings
            return pd.NaT

    elif isinstance(date, cftime.datetime):  # noqa
        for n in range(3):
            try:
                date = pd.Timestamp.fromisoformat((date - pd.Timedelta(n)).isoformat())
            except ValueError:  # We are NOT catching OutOfBoundsDatetime.
                pass
            else:
                break
        else:
            raise ValueError(
                f"Unable to parse cftime date {date}, even when moving back 2 days."
            )
    elif not isinstance(date, pd.Timestamp):
        date = pd.Timestamp(date)  # noqa

    if end_of_period and date_format and not end_date_found:
        if "m" not in date_format:
            date = date + pd.tseries.offsets.YearEnd(1)  # noqa
        elif "d" not in date_format:
            date = date + pd.tseries.offsets.MonthEnd(1)  # noqa

    if output_type == "str":
        return date.strftime(strftime_format)

    return date
