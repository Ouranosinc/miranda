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


def get_ghcn_raw(
    station_ids: list,
    station_type: str,
    out_folder: Path,
    timeout: int = 10,
    update_raw: bool = False,
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

    out_folder.mkdir(parents=True, exist_ok=True)
    errors = []
    for station_id in station_ids:
        if station_type == "daily":
            # url = f"https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/{station_id}.csv"
            url = f"https://noaa-ghcn-pds.s3.amazonaws.com/csv/by_station/{station_id}.csv"

            outfile = out_folder / f"{station_id}.csv"
        # TODO ghcnh not implemented yet
        # elif station_type == "hourly":
        #     url = f"https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/access/by-station/GHCNh_{station_id}_por.psv"
        #     outfile = out_folder / f"GHCNh_{station_id}_por.psv"
        else:
            msg = f"unknown station type : {station_type}"
            raise ValueError(msg)

        if outfile.exists() and not update_raw:
            continue
        try:
            msg = f"Downloading {url}"
            logger.info(msg)
            with requests.get(url, timeout=(5, timeout)) as r:
                r.raise_for_status()
                with Path(outfile.with_suffix(f".tmp{outfile.suffix}")).open("wb") as f:
                    f.write(r.content)
            shutil.move(outfile.with_suffix(f".tmp{outfile.suffix}"), outfile)
        except OSError:
            errors.append(station_id)
            msg = f"Failed to download from URL: {url} .. continuing"
            logger.info(msg)
            continue
    return errors


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
        if project == "ghcnd":
            df = pd.read_csv(station_id)
            df.columns = df.columns.str.lower()
            df.element = df.element.str.lower()
            imask = ~df.q_flag.isin(list(q_flag_dict[project].keys()))
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

                data.append(ds)
        # TODO ghcnh not implemented yet
        # elif project == "ghcnh":
        #     df = pd.read_csv(station_id, delimiter="|")
        #     df.columns = df.columns.str.lower()
        #     # df.element = df.element.str.lower()
        #     # imask = ~df.q_flag.isin(list(q_flag_dict[project].keys()))
        #     df.loc[imask, "q_flag"] = nan
        #     varlist = [k for k in varmeta.keys() if k in df.element.unique()]
        #     if varlist:
        #         df["time"] = pd.to_datetime(df["date"], format="%Y%m%d")
        #         df = df.set_index(["id", "time"])
        #         dslist = []
        #         for var in varlist:
        #             ds1 = df.loc[df.element == var].to_xarray()

        #             ds1 = ds1.rename({"data_value": var, "id": "station"})
        #             drop_vars = [
        #                 v
        #                 for v in ds1.data_vars
        #                 if v not in varlist and v not in ["q_flag"]
        #             ]
        #             ds1 = ds1.drop_vars(drop_vars)
        #             ds1 = ds1.rename(
        #                 {v: f"{var}_{v}" for v in ds1.data_vars if "flag" in v}
        #             )

        #             dslist.append(ds1)
        #         ds = xr.merge(dslist)

        #         del dslist
        #         df_stat = statmeta[statmeta.station_id == station_id.stem]
        #         if len(df_stat) != 1:
        #             raise ValueError(
        #                 f"expected a single station metadata for {station_id.stem}"
        #             )
        #         for cc in [c for c in df_stat.columns if c != "station_id"]:
        #             if cc not in ds.coords:
        #                 ds = ds.assign_coords(
        #                     {
        #                         cc: xr.DataArray(
        #                             df_stat[cc].values, coords=ds.station.coords
        #                         )
        #                     }
        #                 )
        #         for vv in ds.data_vars:
        #             if ds[vv].dtype == "float64":
        #                 ds[vv] = ds[vv].astype("float32")

        #         data.append(ds)
        else:
            msg = f"Unknown project {project}"
            raise ValueError(msg)

    if len(data) == 0:
        return None
    return xr.concat(data, dim="station")


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
        Number of workers to use. Not implemented.

    Raises
    ------
    ValueError
        If the project name is unknown.
    OSError
        If there is an error downloading the data.
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
        )
        if len(errors) == 0:
            break
        else:
            station_ids = errors
            msg = f"Failed to download {len(errors)} stations. Retrying ({try_iter + 1}/{retry})"
            logger.info(msg)
    else:
        msg = f"Failed to download stations after {retry} retries. Giving up."
        logger.error(msg)
        raise RuntimeError(msg)
