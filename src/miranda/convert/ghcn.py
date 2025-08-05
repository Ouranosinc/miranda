"""Module to download and convert GHCN data to Zarr format."""

from __future__ import annotations

import datetime as dt
import logging
import multiprocessing as mp
import os
import shutil
from collections.abc import Generator
from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd
import requests
import xarray as xr
from dask.diagnostics import ProgressBar
from numpy import nan

from miranda.convert._data_corrections import (
    dataset_conversion,
    load_json_data_mappings,
)

logger = logging.getLogger("miranda.convert.ghcn")

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
    }
}

prj_dict = dict(
    ghcnd=dict(freq="daily", filetype=".csv"),
    # ghcnh=dict(freq="hourly", filetype=".psv"), # TODO ghcnh not implemented yet
)


def get_ghcn_raw(
    station_ids: list,
    station_type: str,
    outfolder: Path,
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
    outfolder : Path
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
    if outfolder is None:
        raise ValueError("outfolder must be provided")

    outfolder.mkdir(parents=True, exist_ok=True)
    errors = []
    for station_id in station_ids:
        if station_type == "daily":
            # url = f"https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/{station_id}.csv"
            url = f"https://noaa-ghcn-pds.s3.amazonaws.com/csv/by_station/{station_id}.csv"

            outfile = outfolder / f"{station_id}.csv"
        # TODO ghcnh not implemented yet
        # elif station_type == "hourly":
        #     url = f"https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/access/by-station/GHCNh_{station_id}_por.psv"
        #     outfile = outfolder / f"GHCNh_{station_id}_por.psv"
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


def create_ghcn_xarray(
    infiles: list, varmeta: dict, statmeta: pd.DataFrame, project: str
) -> xr.Dataset | None:
    """
    Create a Zarr dump of DWD climate summary data.

    Parameters
    ----------
    infiles : list
        A list of input files.
    varmeta : dict
        Variable metadata.
    statmeta : pd.DataFrame
        Station metadata.
    project : str
        Project name.

    Returns
    -------
    xr.Dataset, optional
        Dataset.
    """
    data = []
    statmeta
    for station_id in sorted(list(infiles)):
        msg = f"Reading {station_id.name}"
        logger.info(msg)
        if project == "ghcnd":
            df = pd.read_csv(station_id)
            df.columns = df.columns.str.lower()
            df.element = df.element.str.lower()
            imask = ~df.q_flag.isin(list(q_flag_dict[project].keys()))
            df.loc[imask, "q_flag"] = nan
            varlist = [k for k in varmeta.keys() if k in df.element.unique()]
            if varlist:
                df["time"] = pd.to_datetime(df["date"], format="%Y%m%d")
                df = df.set_index(["id", "time"])
                dslist = []
                for var in varlist:
                    ds1 = df.loc[df.element == var].to_xarray()

                    ds1 = ds1.rename({"data_value": var, "id": "station"})
                    drop_vars = [
                        v
                        for v in ds1.data_vars
                        if v not in varlist and v not in ["q_flag"]
                    ]
                    ds1 = ds1.drop_vars(drop_vars)
                    ds1 = ds1.rename(
                        {v: f"{var}_{v}" for v in ds1.data_vars if "flag" in v}
                    )

                    dslist.append(ds1)
                ds = xr.merge(dslist)

                del dslist
                df_stat = statmeta[statmeta.station_id == station_id.stem]
                if len(df_stat) != 1:
                    raise ValueError(
                        f"expected a single station metadata for {station_id.stem}"
                    )
                for cc in [
                    c
                    for c in df_stat.columns
                    if c not in ["station_id", "geometry", "index_right"]
                ]:
                    if cc not in ds.coords:
                        ds = ds.assign_coords(
                            {
                                cc: xr.DataArray(
                                    [df_stat[cc].values[0]], coords=ds.station.coords
                                )
                            }
                        )
                for vv in ds.data_vars:
                    if ds[vv].dtype == "float64":
                        ds[vv] = ds[vv].astype("float32")

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
    n_workers: int | None = None,
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
    n_workers : int, optional
        Number of workers to use. Not implemented.
    """
    station_df = _get_ghcn_stations(
        project=project, lon_bnds=lon_bnds, lat_bnds=lat_bnds
    )
    if update_raw and working_folder.joinpath("raw").exists():
        shutil.rmtree(working_folder.joinpath("raw"))
    working_folder.mkdir(parents=True, exist_ok=True)
    working_folder.joinpath("raw").mkdir(exist_ok=True)

    # request = NoaaGhcnRequest(parameters=(prj_dict[project], "data"), start_date=start_date, end_date=end_date)

    station_ids = station_df["station_id"].tolist()

    ntry = 5
    while ntry > 0:
        errors = get_ghcn_raw(
            station_ids=station_ids,
            station_type=prj_dict[project]["freq"],
            outfolder=working_folder.joinpath("raw"),
            update_raw=update_raw,
            timeout=timeout,
        )
        if len(errors) == 0:
            break
        else:
            station_ids = errors
            ntry -= 1
            msg = f"Failed to download {len(errors)} stations. Retrying ntry={ntry}"
            logger.info(msg)


def _get_ghcn_stations(
    project: str,
    lon_bnds: list[float] | None = None,
    lat_bnds: list[float] | None = None,
) -> pd.DataFrame:
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


def convert_ghcn_bychunks(
    project: str,
    working_folder: str | os.PathLike[str] | None = None,
    cfvariable_list: list | None = None,
    start_year: int | None = None,
    end_year: int | None = None,
    lon_bnds: list[float] | None = None,
    lat_bnds: list[float] | None = None,
    n_workers: int = 4,
    nstations: int = 100,
    update_from_raw: bool = False,
) -> None:
    """
    Convert GHCN data to Zarr format.

    Requires GIS libraries (geopandas).

    Parameters
    ----------
    project : str
        Project name.
    working_folder : str or os.PathLike[str], optional
        The working folder. The default (None) is to use the current working directory.
    cfvariable_list : list, optional
        List of CF variable names. Optional.
    start_year : int, optional
        Start year. Optional.
    end_year : int, optional
        End year. Optional.
    lon_bnds : list of float, optional
        Longitude boundaries.
    lat_bnds : list of float, optional
        Latitude boundaries.
    n_workers : int
        Number of workers to use. Default is 4.
    nstations : int
        Number of stations to process. Default is 100.
    update_from_raw : bool
        Whether to update from raw data.
    """
    try:
        import geopandas as gpd
        from shapely.geometry import box
    except ImportError:
        msg = "GNCN conversion requires the GIS libraries. Install them with `$ pip install miranda[gis]`."
        logger.error(msg)
        raise

    var_attrs = load_json_data_mappings(project=project)["variables"]
    if cfvariable_list:
        var_attrs = {
            v: var_attrs[v]
            for v in var_attrs
            if var_attrs[v]["_cf_variable_name"] in cfvariable_list
        }
    freq_dict = dict(h="hr", d="day")
    if project == "ghcnd":

        readme_url = "https://noaa-ghcn-pds.s3.amazonaws.com/readme.txt"

        outchunks = dict(time=(365 * 4) + 1, station=nstations)
    # TODO ghcnh not implemented yet
    # elif project == "ghcnh":
    #     readme_url = "https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/doc/ghcnh_DOCUMENTATION.pdf"
    #     outchunks = dict(time=(365 * 4) + 1, station=nstations)
    # logger.info("ghcnh not implemented yet")
    # exit()
    else:
        msg = f"Unknown project {project}"
        raise ValueError(msg)
    station_df = _get_ghcn_stations(
        project=project, lon_bnds=lon_bnds, lat_bnds=lat_bnds
    )
    tz_file = Path(__file__).parent.joinpath(
        "data/timezones-with-oceans-now.shapefile.zip"
    )

    tz = gpd.read_file(tz_file).to_crs(epsg=4326)
    # clip to bbox for faster sjoin
    # Create a custom polygon
    polygon = box(
        lon_bnds[0] - 0.1, lat_bnds[0] - 0.1, lon_bnds[-1] + 0.1, lat_bnds[-1] + 0.1
    )
    poly_clip = gpd.GeoDataFrame([1], geometry=[polygon], crs=tz.crs)
    tz = tz.clip(poly_clip)
    station_df = gpd.GeoDataFrame(
        station_df,
        geometry=gpd.points_from_xy(station_df.lon, station_df.lat),
        crs=tz.crs,
    ).sjoin(tz, how="left")
    station_df = station_df.rename(columns={"tzid": "timezone"})
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()
    working_folder.mkdir(parents=True, exist_ok=True)
    working_folder.joinpath("zarr").mkdir(exist_ok=True)

    if end_year is not None:
        end_date = dt.datetime(end_year, 12, 31, 23, 59, 59, tzinfo=ZoneInfo("UTC"))
    else:
        end_date = dt.datetime.now().astimezone(ZoneInfo("UTC"))
    start_date = dt.datetime(start_year, 1, 1, 0, 0, 0, tzinfo=ZoneInfo("UTC"))

    if update_from_raw:
        for folder in working_folder.joinpath("zarr").iterdir():
            msg = f"Deleting {folder}"
            logger.info(msg)
            shutil.rmtree(folder)

    def _chunk_list(lst: list, n: int) -> Generator[list]:
        """Split list into chunks of size n."""
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    treated = []
    file_list = sorted(
        list(working_folder.joinpath("raw").rglob(f"*{prj_dict[project]['filetype']}"))
    )
    jobs = []
    for ii, ss in enumerate(_chunk_list(file_list, nstations)):
        if ii not in treated:
            var_attrs_new = {}
            for vv, meta in var_attrs.items():
                cf_var = var_attrs[vv]["_cf_variable_name"]
                outzarr = working_folder.joinpath(
                    "zarr", cf_var, f"{project}_{ii}.zarr"
                )
                if not outzarr.exists() or update_from_raw:
                    var_attrs_new[vv] = meta

            if var_attrs_new:
                dsall_vars = create_ghcn_xarray(
                    infiles=ss,
                    varmeta=var_attrs_new,
                    statmeta=station_df,
                    project=project,
                )
                if dsall_vars is None:
                    continue
                dsall_vars = dsall_vars.sel(
                    time=slice(str(start_date.year), str(end_date.year))
                )
                for kk, vv in var_attrs_new.items():
                    cf_var = var_attrs[kk]["_cf_variable_name"]
                    outzarr = working_folder.joinpath(
                        "zarr", cf_var, f"{project}_{ii}.zarr"
                    )
                    if kk not in dsall_vars.data_vars:
                        continue
                    dsout = dsall_vars.drop_vars(
                        [v for v in dsall_vars.data_vars if not v.startswith(kk)]
                    )
                    allnull_stat = dsout[kk].isnull().sum(dim="time") == len(dsout.time)
                    dsout = dsout.sel(station=~allnull_stat)
                    dsout = make_monotonous_time(dsout, freq=prj_dict[project]["freq"])

                    ds_corr = dataset_conversion(
                        dsout,
                        project=project,
                        add_version_hashes=False,
                        overwrite=update_from_raw,
                    )
                    ds_corr = ds_corr.rename({f"{kk}_q_flag": f"{cf_var}_q_flag"})
                    for vv in ds_corr.data_vars:
                        if ds_corr[vv].dtype == "float64":
                            ds_corr[vv] = ds_corr[vv].astype("float32")

                    desc_str = "; ".join(
                        [f"{k}:{v}" for k, v in q_flag_dict[project].items()]
                    )
                    desc_str = f"{desc_str}. See the readme file for information of quality flag (QFLAG1) codes : {readme_url}"
                    attrs = {
                        "flag_values": [c for c in q_flag_dict[project].keys()],
                        "flag_meanings": [c for c in q_flag_dict[project].values()],
                        "standard_name": f"{ds_corr[cf_var].attrs['standard_name']}_quality_flag",
                        "long_name": f"Quality flag for {cf_var}",
                        "description": desc_str,
                    }

                    ds_corr[f"{cf_var}_q_flag"].attrs = attrs

                    jobs.append((ds_corr, outzarr, outchunks))
                    if len(jobs) >= n_workers:
                        pool = mp.Pool(n_workers)
                        pool.starmap(write_zarr, jobs)
                        pool.close()
                        pool.join()
                        jobs = []

            if len(jobs) > 0:
                pool = mp.Pool(n_workers)
                pool.starmap(write_zarr, jobs)
                pool.close()
                pool.join()
                jobs = []
            treated.append(ii)


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
    outzarr: Path,
    chunks: dict,
    overwrite: bool = False,
) -> None:
    """
    Write Zarr file.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset.
    outzarr : Path
        Output Zarr file.
    chunks : dict
        Chunk sizes.
    overwrite : bool
        Whether to overwrite. Default is False.
    """
    if not outzarr.exists() or overwrite:
        with ProgressBar():
            ds.chunk(chunks).to_zarr(outzarr.with_suffix(".tmp.zarr"), mode="w")
        shutil.move(outzarr.with_suffix(".tmp.zarr"), outzarr)
