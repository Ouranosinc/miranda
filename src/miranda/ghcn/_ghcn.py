from __future__ import annotations
import logging
import os
import shutil
from pathlib import Path

import pandas as pd
import requests
import xarray as xr
from numpy import nan

from miranda.convert.utils import (
    _add_coords_to_dataset,
    get_station_meta,
    prj_dict,
    q_flag_dict,
)


all = [
    "create_ghcn_xarray",
    "download_ghcn",
    "get_ghcn_raw",
]

logger = logging.getLogger("miranda.ghcn")


def _process_ghcnd(station_id: Path, variable_meta: dict, station_meta: pd.DataFrame, logger: logging.Logger) -> xr.Dataset | None:
    """
    Process a single GHCN-Daily (ghcnd) station file into an xarray.Dataset.

    Parameters
    ----------
    station_id : Path
        Path to the input CSV file for the station.
    variable_meta : dict
        Variable metadata dictionary.
    station_meta : pd.DataFrame
        DataFrame containing station metadata.
    logger : logging.Logger
        Logger for logging messages.

    Returns
    -------
    xr.Dataset or None
        The processed dataset, or None if no variables found.
    """
    df = pd.read_csv(station_id)
    df.columns = df.columns.str.lower()
    df.element = df.element.str.lower()
    imask = ~df.q_flag.isin(list(q_flag_dict["ghcnd"].keys()))
    df.loc[imask, "q_flag"] = nan
    varlist = [k for k in variable_meta.keys() if k in df.element.unique()]
    if varlist:
        df["time"] = pd.to_datetime(df["date"], format="%Y%m%d")
        df = df.set_index(["id", "time"])
        dslist = []
        for var in varlist:
            ds1 = df.loc[df.element == var].to_xarray()
            ds1 = ds1.rename({"data_value": var, "id": "station"})
            drop_vars = [v for v in ds1.data_vars if v not in varlist and v not in ["q_flag"]]
            ds1 = ds1.drop_vars(drop_vars)
            ds1 = ds1.rename({v: f"{var}_{v}" for v in ds1.data_vars if "flag" in v})
            dslist.append(ds1)
        ds = xr.merge(dslist)
        del dslist
        df_stat = station_meta[station_meta.station_id == station_id.stem]
        if len(df_stat) != 1:
            raise ValueError(f"expected a single station metadata for {station_id.stem}")
        ds = _add_coords_to_dataset(ds, df_stat, float_flag=False)
        return ds
    return None


def _process_ghcnh(station_id: Path, variable_meta: dict, station_meta: pd.DataFrame, logger: logging.Logger) -> xr.Dataset | None:
    """
    Process a single GHCN-Hourly (ghcnh) station file into an xarray.Dataset.

    Parameters
    ----------
    station_id : Path
        Path to the input PSV file for the station.
    variable_meta : dict
        Variable metadata dictionary.
    station_meta : pd.DataFrame
        DataFrame containing station metadata.
    logger : logging.Logger
        Logger for logging messages.

    Returns
    -------
    xr.Dataset or None
        The processed dataset, or None if no variables found.
    """
    df = pd.read_csv(station_id, delimiter="|", low_memory=False)
    df.columns = df.columns.str.lower()
    varlist = [k for k in variable_meta.keys() if k in df.columns]
    flaglist = [f"{k}_quality_code" for k in varlist]
    coordlist = [
        "station_id",
        "station_name",
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "latitude",
        "longitude",
        "elevation",
    ]
    drop_cols = [c for c in df.columns if c not in varlist and c not in flaglist and c not in coordlist]
    df = df.drop(columns=drop_cols)
    if varlist:
        for col in ["year", "month", "day", "hour"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            df = df.dropna(subset=[col])
        df["time"] = pd.to_datetime(df[["year", "month", "day", "hour"]])
        df = df.set_index(["station_id", "time"])
        df = df.iloc[~df.index.duplicated()]
        dslist = []
        for var in varlist:
            ds1 = df[[var, f"{var}_quality_code"]].to_xarray()
            ds1 = ds1.rename({"station_id": "station", f"{var}_quality_code": f"{var}_flag"})
            if ds1[f"{var}_flag"].dtype == "float":
                ds1[f"{var}_flag"] = ds1[f"{var}_flag"].round().astype(str)
                ds1[f"{var}_flag"] = ds1[f"{var}_flag"].where(ds1[f"{var}_flag"] != "nan", "")
            else:
                df1 = df[[var, f"{var}_quality_code"]].copy()
                df1["num_str"] = pd.to_numeric(df1[f"{var}_quality_code"], errors="coerce").round().astype(str)
                df1.loc[
                    df1[f"{var}_quality_code"].apply(lambda x: isinstance(x, str)),
                    "num_str",
                ] = df1[df1[f"{var}_quality_code"].apply(lambda x: isinstance(x, str))][f"{var}_quality_code"].astype(str)
                df1 = df1.drop(columns=[f"{var}_quality_code"])
                df1 = df1.rename(columns={"num_str": f"{var}_quality_code"})
                ds1 = df1.to_xarray()
                ds1 = ds1.rename({"station_id": "station", f"{var}_quality_code": f"{var}_flag"})
                ds1[f"{var}_flag"] = ds1[f"{var}_flag"].where(ds1[f"{var}_flag"] != "nan", "")
            dslist.append(ds1)
        ds = xr.merge(dslist)
        del dslist
        df_stat = station_meta[station_meta.station_id == station_id.stem.split("_", 1)[1].split("_por")[0]]
        if len(df_stat) != 1:
            raise ValueError(f"expected a single station metadata for {station_id.stem}")
        for cc in [c for c in df_stat.columns if c not in ["station_id", "geometry", "index_right"]]:
            if cc not in ds.coords:
                ds = ds.assign_coords({cc: xr.DataArray(df_stat[cc].values, coords=ds.station.coords)})
        for vv in ds.data_vars:
            if ds[vv].dtype == "float64":
                ds[vv] = ds[vv].astype("float32")
        return ds
    return None


def create_ghcn_xarray(in_files: list, variable_meta: dict, station_meta: pd.DataFrame, project: str) -> xr.Dataset | None:
    """
    Create a Zarr dump of DWD climate summary data.

    Parameters
    ----------
    in_files : list
        A list of input files.
    variable_meta : dict
        Variable metadata.
    station_meta : pd.DataFrame
        Station metadata.
    project : str
        Project name.

    Returns
    -------
    xr.Dataset, optional
        Dataset.
    """
    data = []
    for station_id in sorted(list(in_files)):
        msg = f"Reading {station_id.name}"
        logger.info(msg)
        if project in prj_dict:
            try:
                if project == "ghcnd":
                    ds = _process_ghcnd(station_id, variable_meta, station_meta, logger)
                elif project == "ghcnh":
                    ds = _process_ghcnh(station_id, variable_meta, station_meta, logger)
                else:
                    msg = f"Unknown project {project}"
                    raise ValueError(msg)
                if ds is not None:
                    data.append(ds)
            except Exception as e:
                msg = f"Failed to read data for {station_id.name} : {e} ... continuing"
                logger.warning(msg)
                continue
        else:
            msg = f"Unknown project {project}"
            raise ValueError(msg)

    if len(data) == 0:
        return None
    return xr.concat(data, dim="station", join="outer")


def get_ghcn_raw(
    station_ids: list,
    station_type: str,
    out_folder: Path,
    timeout: int = 10,
    update_raw: bool = False,
    n_workers: int | None = None,
) -> list[str]:
    """
    Download raw GHCN data.

    Parameters
    ----------
    station_ids : list[str]
        List of station IDs.
    station_type : str
        Station type.
    out_folder : Path
        Output folder.
    timeout : int
        Request timeout in seconds. Default is 10.
    update_raw : bool
        Whether to update raw data.
    n_workers : int, optional
        Number of parallel workers to use. If None or 1, no parallelism is used

    Returns
    -------
    list of str
        List of station IDs that failed to download.
    """
    if station_ids is None:
        raise ValueError("station_ids must be provided")
    if station_type is None:
        raise ValueError("stationtype must be provided")
    if out_folder is None:
        raise ValueError("outfolder must be provided")

    import concurrent.futures

    out_folder.mkdir(parents=True, exist_ok=True)
    errors = []

    def download_one_station(station_id):
        if station_type == "daily":
            url = f"https://noaa-ghcn-pds.s3.amazonaws.com/csv/by_station/{station_id}.csv"
            outfile = out_folder / f"{station_id}.csv"
        elif station_type == "hourly":
            url = f"https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/access/by-station/GHCNh_{station_id}_por.psv"
            outfile = out_folder / f"GHCNh_{station_id}_por.psv"
        else:
            msg = f"unknown station type : {station_type}"
            raise ValueError(msg)

        if outfile.exists() and not update_raw:
            return None
        try:
            msg = f"Downloading {url}"
            logger.info(msg)
            with requests.get(url, timeout=(5, timeout)) as r:
                r.raise_for_status()
                with Path(outfile.with_suffix(f".tmp{outfile.suffix}")).open("wb") as f:
                    f.write(r.content)
            shutil.move(outfile.with_suffix(f".tmp{outfile.suffix}"), outfile)
            return None
        except OSError:
            msg = f"Failed to download from URL: {url} .. continuing"
            logger.info(msg)
            return station_id

    if n_workers is not None and n_workers > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
            results = list(executor.map(download_one_station, station_ids))
        errors = [sid for sid in results if sid is not None]
    else:
        for station_id in station_ids:
            err = download_one_station(station_id)
            if err is not None:
                errors.append(err)
    return errors


def download_ghcn(
    project: str,
    working_folder: str | os.PathLike[str] | None = None,
    lon_bnds: list[float] | None = None,
    lat_bnds: list[float] | None = None,
    update_raw: bool = False,
    timeout: int | None = None,
    retry: int = 5,
    n_workers: int | None = None,  # FIXME: Not implemented yet
) -> None:
    """
    Download GHCN data.

    Parameters
    ----------
    project : str
        Project name.
    working_folder : str or os.PathLink[str], optional
        Temporary files folder.
    lon_bnds : list of float, optional
        Longitude boundaries.
    lat_bnds : list of float, optional
        Latitude boundaries.
    update_raw : bool
        Whether to update the raw files or not.
    timeout : int, optional
        Request timeout in seconds.
    retry : int
        Number of retries.
    n_workers : int, optional
        Number of parallel workers to use. If None or 1, no parallelism is used

    Raises
    ------
    ValueError
        If the project name is unknown.
    """
    station_df = get_station_meta(project=project, lon_bnds=lon_bnds, lat_bnds=lat_bnds)
    if update_raw and working_folder.joinpath("raw").exists():
        shutil.rmtree(working_folder.joinpath("raw"))
    working_folder.mkdir(parents=True, exist_ok=True)

    out_folder = working_folder.joinpath("raw")
    out_folder.mkdir(parents=True, exist_ok=True)

    # request = NoaaGhcnRequest(parameters=(prj_dict[project], "data"), start_date=start_date, end_date=end_date)

    station_ids = station_df["station_id"].tolist()

    for try_iter in range(retry):
        errors = get_ghcn_raw(
            station_ids=station_ids,
            station_type=prj_dict[project]["freq"],
            out_folder=out_folder,
            update_raw=update_raw,
            timeout=timeout,
            n_workers=n_workers,
        )
        if len(errors) == 0:
            break
        else:
            station_ids = errors
            msg = f"Failed to download {len(errors)} stations. Retrying ({try_iter + 1}/{retry})"
            logger.info(msg)
    else:
        msg = f"Failed to download stations {errors} after {retry} retries....skipping."
        logger.error(msg)
