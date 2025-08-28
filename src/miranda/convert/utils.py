"""Conversion Utilities submodule."""

from __future__ import annotations
import hashlib
import logging
import os
import re
import shutil
from pathlib import Path
from typing import Any

import cftime
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar
from pandas._libs import NaTType  # noqa


logger = logging.getLogger("miranda.convert.utils")

__all__ = [
    "_add_coords_to_dataset",
    "date_parser",
    "find_version_hash",
    "get_station_meta",
]

prj_dict = dict(
    ghcnd=dict(freq="daily", filetype=".csv"),
    canhomt_dly=dict(freq="daily", filetype=".csv"),
    # ghcnh=dict(freq="hourly", filetype=".psv"), # TODO ghcnh not implemented yet
)

q_flag_dict = {
    "ghcnd": {
        "D": "failed duplicate check",
        "G": "failed gap check",
        "I": "failed internal consistency check",
        "K": "failed streak/frequent-value check",
        "L": "failed check on length of multiday period",
        "M": "failed megaconsistency check",
        "N": "failed naught check",
        "O": "failed climatological outlier check",
        "R": "failed lagged range check",
        "S": "failed spatial consistency check",
        "T": "failed temporal consistency check",
        "W": "temperature too warm for snow",
        "X": "failed bounds check",
        "Z": "flagged as a result of an official Datzilla investigation",
    },
    "canhomt_dly": {
        "I": "infilled data",
    },
}


def _add_coords_to_dataset(ds: xr.Dataset, df_stat: pd.DataFrame, float_flag=True) -> xr.Dataset:
    """
    Add coordinates to the dataset from the station metadata.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset.
    df_stat : pd.DataFrame
        Station metadata.
    float_flag : bool, optional
        Whether to convert data variables to float32. Default is True.

    Returns
    -------
    xr.Dataset
        Dataset with added coordinates.
    """
    for cc in [c for c in df_stat.columns if c not in ["station_id", "geometry", "index_right"]]:
        if cc not in ds.coords:
            ds = ds.assign_coords({cc: xr.DataArray([df_stat[cc].values[0]], coords=ds.station.coords)})
    if float_flag:
        for vv in ds.data_vars:
            if ds[vv].dtype == "float64":
                ds[vv] = ds[vv].astype("float32")
    return ds


def find_version_hash(file: str | os.PathLike[str]) -> dict[str, Any]:
    """
    Check for an existing version hash file and, if one cannot be found, generate one from file.

    Parameters
    ----------
    file : str or os.PathLike
        The file to check.

    Returns
    -------
    dict
        The version and hash.
    """

    def _get_hash(f: str) -> str:
        """
        Calculate the sha256sum of a file.

        Parameters
        ----------
        f : str or os.PathLike
            The file to hash.

        Returns
        -------
        str
            The hash.
        """
        hash_sha256_writer = hashlib.sha256()
        with Path(f).open("rb") as f_opened:
            hash_sha256_writer.update(f_opened.read())
        sha256sum = hash_sha256_writer.hexdigest()
        _msg = f"Calculated sha256sum (starting: {sha256sum[:6]})"
        logger.info(_msg)
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
                    msg = "Unable to read version hash file. Calculating sha256sum."
                    logger.error(msg)
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
) -> str | pd.Timestamp | NaTType:
    """
    Parse datetime objects from a string representation of a date or both a start and end date.

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
    Adapted from code written by Gabriel Rondeau-Genesse (@RondeauG).
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

    def _parse_date(d: str, fmts: list[str]) -> tuple[pd.Timestamp, str]:
        """
        Parse the date.

        Parameters
        ----------
        d : str
            The date string.
        fmts : list
            The list of formats to try.

        Returns
        -------
        pd.Timestamp
            The parsed date.
        """
        for fmt in fmts:
            try:
                s = pd.to_datetime(d, format=fmt)
                match = fmt
                break
            except ValueError:  # noqa: S110
                pass
        else:
            raise ValueError(f"Can't parse date {d} with supported formats: [{', '.join(fmts)}].")
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
            except (  # noqa: PERF203,S110  # We are NOT catching OutOfBoundsDatetime.
                ValueError
            ):
                pass
            else:
                break
        else:
            raise ValueError(f"Unable to parse cftime date {date}, even when moving back 2 days.")
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


def get_station_meta(
    project: str,
    lon_bnds: list[float] | None = None,
    lat_bnds: list[float] | None = None,
) -> pd.DataFrame:
    """
    Get GHCN or CanHomT station metadata.

    Parameters
    ----------
    project : str
        Project name.
    lon_bnds : list of float, optional
        Longitude boundaries.
    lat_bnds : list of float, optional
        Latitude boundaries.

    Returns
    -------
    pd.DataFrame
        Station metadata.
    """
    if project == "ghcnd":
        station_df = _get_ghcn_stations(project=project)
    elif project == "canhomt_dly":
        station_df = _get_canhomt_stations(project=project)

    # TODO ghcnh not implemented yet
    # elif project == "ghcnh":
    #     project
    #     station_url = "https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/doc/ghcnh-station-list.txt"
    #     dtypes = {
    #         "station_id": str,
    #         "lat": float,
    #         "lon": float,
    #         "elevation": float,
    #         "state": str,
    #         "station_name": str,
    #         "gsn_flag": str,
    #         "hcn_flag": str,
    #         "wmo_id": str,
    #     }
    #     try:
    #         station_df = pd.read_fwf(
    #             station_url,
    #             widths=[11, 9, 10, 7, 3, 31, 4, 4, 6],
    #             names=[d for d in dtypes.keys()],
    #             converters=dtypes,
    #         )
    #     except ValueError:
    #         statfile = Path(__file__).parent.joinpath("data/ghcnh-station-list.txt")
    #         station_df = pd.read_fwf(
    #             statfile,
    #             widths=[11, 9, 10, 7, 3, 31, 4, 4, 6],
    #             names=[d for d in dtypes.keys()],
    #             converters=dtypes,
    #         )

    # logger.info("ghcnh not implemented yet")
    # exit()
    else:
        raise ValueError(f"unknown project values {project}")
    if lon_bnds:  # and lat_bnds:
        bbx_mask = station_df["lon"].between(lon_bnds[0], lon_bnds[1])
        station_df = station_df[bbx_mask]
    if lat_bnds:
        bbx_mask = station_df["lat"].between(lat_bnds[0], lat_bnds[1])
        station_df = station_df[bbx_mask]

    return station_df


# FIXME: Needs to be moved to different module in a future version
def _get_ghcn_stations(project: str) -> pd.DataFrame:
    """
    Get GHCN station metadata.

    Parameters
    ----------
    project : str
        Project name.

    Returns
    -------
    pd.DataFrame
        Station metadata.
    """
    if project == "ghcnd":
        station_url = "https://noaa-ghcn-pds.s3.amazonaws.com/ghcnd-stations.txt"
        dtypes = {
            "station_id": str,
            "lat": float,
            "lon": float,
            "elevation": float,
            "state": str,
            "station_name": str,
            "gsn_flag": str,
            "hcn_flag": str,
            "wmo_id": str,
        }
        try:
            station_df = pd.read_fwf(
                station_url,
                widths=[11, 9, 10, 7, 3, 31, 4, 4, 6],
                names=[d for d in dtypes.keys()],
                converters=dtypes,
            )
        except ValueError:
            statfile = Path(__file__).parent.joinpath("data/ghcnd-stations.txt")
            station_df = pd.read_fwf(
                statfile,
                widths=[11, 9, 10, 7, 3, 31, 4, 4, 6],
                names=[d for d in dtypes.keys()],
                converters=dtypes,
            )
    else:
        raise ValueError(f"unknown project values {project}")
    return station_df


# FIXME: Needs to be moved to different module in a future version
def _get_canhomt_stations(project: str) -> pd.DataFrame:
    """
    Get CanHomT station metadata.

    Parameters
    ----------
    project : str
        Project name.

    Returns
    -------
    pd.DataFrame
        Station metadata.
    """
    if project == "canhomt_dly":
        station_url = "https://crd-data-donnees-rdc.ec.gc.ca/CDAS/products/CanHomTV4/Temp_Stations_Gen4_2024_monthly.csv"

        try:
            station_df = pd.read_csv(
                station_url,
            )
        except ValueError:
            statfile = Path(__file__).parent.joinpath("data/eccc-canhomt_Temp_Stations_Gen4_2024_monthly.zip")
            station_df = pd.read_csv(
                statfile,
            )
        station_df = station_df.rename(columns={p: p.lower() for p in station_df.columns})
        rename = {"name": "station_name", "id": "station_id", "ele": "elevation"}
        station_df = station_df.rename(columns=rename)
    else:
        raise ValueError(f"Unknown project value: {project}")

    return station_df


def make_monotonous_time(ds: xr.Dataset, freq: str):
    """
    Make time monotonous.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset.
    freq : str
        Frequency.

    Returns
    -------
    xr.Dataset
        Dataset with monotonous time.
    """
    if freq == "daily":
        time1 = pd.date_range(ds.time[0].values, ds.time[-1].values, freq="D")
    elif freq == "hourly":
        time1 = pd.date_range(ds.time[0].values, ds.time[-1].values, freq="H")
    else:
        raise ValueError(f"Unknown frequency {freq}")
    dsnew = xr.Dataset(coords=dict(time=time1, station=ds.station))
    for vv in ds.data_vars:
        dsnew[vv] = ds[vv]
    return dsnew


def write_zarr(
    ds: xr.Dataset,
    out_zarr: Path,
    chunks: dict,
    overwrite: bool = False,
) -> None:
    """
    Write Zarr file.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset.
    out_zarr : Path
        Output Zarr file.
    chunks : dict
        Chunk sizes.
    overwrite : bool
        Whether to overwrite. Default is False.
    """
    if not out_zarr.exists() or overwrite:
        with ProgressBar():
            ds.chunk(chunks).to_zarr(out_zarr.with_suffix(".tmp.zarr"), mode="w")
        shutil.move(out_zarr.with_suffix(".tmp.zarr"), out_zarr)
