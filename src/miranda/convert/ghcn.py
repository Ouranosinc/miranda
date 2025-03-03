from __future__ import annotations

import datetime as dt
import logging.config
import multiprocessing as mp
import os
import shutil
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import pandas as pd
import requests
import xarray as xr
from dask.diagnostics import ProgressBar
from numpy import unique

from miranda.convert._data_corrections import (
    dataset_conversion,
    load_json_data_mappings,
)
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


def chunk_list(lst, n):
    """Split list into chunks of size n."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def get_ghcn_raw(
    station_ids: list = None,
    stationtype: str = None,
    outfolder: Path = None,
    timeout: int = None,
    update_raw: bool = False,
) -> None:
    """Download raw GHCN data."""
    if station_ids is None:
        raise ValueError("station_ids must be provided")
    if stationtype is None:
        raise ValueError("stationtype must be provided")
    if outfolder is None:
        raise ValueError("outfolder must be provided")
    if timeout is None:
        timeout = 10

    outfolder.mkdir(parents=True, exist_ok=True)
    errors = []
    for station_id in station_ids:
        if stationtype == "daily":
            # url = f"https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/{station_id}.csv"
            url = f"https://noaa-ghcn-pds.s3.amazonaws.com/csv/by_station/{station_id}.csv"

            outfile = outfolder / f"{station_id}.csv"
            if outfile.exists() and not update_raw:
                continue
            try:
                logging.info(f"downloading {url}")
                with requests.get(url, timeout=(5, timeout)) as r:
                    r.raise_for_status()
                    with open(outfile.with_suffix(".tmp.csv"), "wb") as f:
                        f.write(r.content)
                shutil.move(outfile.with_suffix(".tmp.csv"), outfile)
            except:
                errors.append(station_id)
                logging.info(f"Failed to download {url} .. continuing")
                continue
    return errors


def create_ghcn_xarray(infolder: Path, varmeta: dict, statmeta: pd.DataFrame) -> None:
    """Create a Zarr dump of DWD climate summary data."""
    data = []
    statmeta
    for station_id in sorted(list(infolder.glob("*.csv"))):
        logging.info(f"reading {station_id.name}")
        df = pd.read_csv(station_id)
        df.columns = df.columns.str.lower()
        df.element = df.element.str.lower()
        inull = df["q_flag"].isnull()
        df["q_flag"] = df["q_flag"].astype(str)
        df.loc[inull, "q_flag"] = ""
        df.loc[df["q_flag"] == "nan", "q_flag"] = ""
        varlist = [k for k in varmeta.keys() if k in df.element.unique()]
        if varlist:

            # cntryid = [x[1].split(' ')[-1] for x in df['name'].str.split(', ')]
            # stateid = [x[1].split(' ')[0] for x in df['name'].str.split(', ')]
            # stat_name = [x[0] for x in df['name'].str.split(', ')]
            # df['stateid'] = stateid
            # df['cntryid'] = cntryid
            # df['station_name'] = stat_name

            df["time"] = pd.to_datetime(df["date"], format="%Y%m%d")
            df = df.set_index(["id", "time"])
            dslist = []
            for var in varlist:
                ds1 = df.loc[df.element == var].to_xarray()
                ds1 = ds1.rename({"data_value": var, "id": "station"})
                drop_vars = [
                    v for v in ds1.data_vars if v not in varlist and v not in ["q_flag"]
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
            for cc in [c for c in df_stat.columns if c != "station_id"]:
                if cc not in ds.coords:
                    ds = ds.assign_coords(
                        {cc: xr.DataArray(df_stat[cc].values, coords=ds.station.coords)}
                    )
            for vv in ds.data_vars:
                if ds[vv].dtype == "float64":
                    ds[vv] = ds[vv].astype("float32")
                if "flag" in vv:
                    ds[vv] = ds[vv].astype("str")

            data.append(ds)
    if len(data) == 0:
        return None
    return xr.concat(data, dim="station")


def convert_ghcn_bychunks(
    project: str,
    working_folder: str | os.PathLike[str] | None = None,
    overwrite: bool = False,
    cfvariable_list: list | None = None,
    lon_bnds: list[float] | None = None,
    lat_bnds: list[float] | None = None,
    start_year: int | None = None,
    end_year: int | None = None,
    n_workers: int = 4,
    nstations: int = 100,
    update_raw: bool = False,
    delete_raw: bool = False,
) -> None:

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
        ctry_url = "https://noaa-ghcn-pds.s3.amazonaws.com/ghcnd-countries.txt"
        try:
            ctry_df = pd.read_fwf(ctry_url, widths=[2, 50])
        except:
            ctry_file = Path(__file__).parent.joinpath("data/ghcnd-countries.txt")
            ctry_df = pd.read_fwf(ctry_file, widths=[2, 50])

        ctry_df.columns = ["code", "name"]
        ctry_df = ctry_df.set_index("code")
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
        except:
            statfile = Path(__file__).parent.joinpath("data/ghcnd-stations.txt")
            station_df = pd.read_fwf(
                statfile,
                widths=[11, 9, 10, 7, 3, 31, 4, 4, 6],
                names=[d for d in dtypes.keys()],
                converters=dtypes,
            )

    elif project == "ghcnh":
        logging.info("ghcnh not implemented yet")
        exit()

    prj_dict = dict(ghcnd="daily", ghcnh="hourly")
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()
    working_folder.mkdir(parents=True, exist_ok=True)
    working_folder.joinpath("raw").mkdir(exist_ok=True)
    working_folder.joinpath("zarr").mkdir(exist_ok=True)

    bbox = None
    if lon_bnds and lat_bnds:
        bbx_mask = station_df["lat"].between(lat_bnds[0], lat_bnds[1]) & station_df[
            "lon"
        ].between(lon_bnds[0], lon_bnds[1])
    else:
        raise ValueError("latitude and longitude bounds must be provided")
    if end_year is not None:
        end_date = dt.datetime(end_year, 12, 31, 23, 59, 59, tzinfo=ZoneInfo("UTC"))
    else:
        end_date = dt.datetime.now().astimezone(ZoneInfo("UTC"))
    start_date = dt.datetime(start_year, 1, 1, 0, 0, 0, tzinfo=ZoneInfo("UTC"))

    # request = NoaaGhcnRequest(parameters=(prj_dict[project], "data"), start_date=start_date, end_date=end_date)
    station_df = station_df[bbx_mask]
    station_ids = station_df["station_id"].tolist()
    station_list = sorted(list(chunk_list(station_ids, nstations)))
    treated = []

    if update_raw:
        for folder in working_folder.joinpath("raw").iterdir():
            logging.info(f"deleting {folder}")
            shutil.rmtree(folder)
        for folder in working_folder.joinpath("zarr").iterdir():
            logging.info(f"deleting {folder}")
            shutil.rmtree(folder)

    for ii, ss in enumerate(station_list):

        ntry = 5
        while ntry > 0:
            errors = get_ghcn_raw(
                station_ids=ss,
                stationtype=prj_dict[project],
                outfolder=working_folder.joinpath("raw", str(ii)),
                update_raw=update_raw,
            )
            if len(errors) == 0:
                break
            else:
                ss = errors
                ntry -= 1
                logging.info(f"Retrying ntry={ntry}")

        infolders = [
            i
            for i in working_folder.joinpath("raw").iterdir()
            if i.is_dir() and i not in treated
        ]
        if len(infolders) >= n_workers or ii + 1 == len(station_list):

            jobs = []
            for infolder in infolders:

                var_attrs_new = {}
                for vv, meta in var_attrs.items():
                    cf_var = var_attrs[vv]["_cf_variable_name"]
                    outzarr = working_folder.joinpath(
                        "zarr", cf_var, f"{project}_{infolder.name}.zarr"
                    )
                    if not outzarr.exists() or overwrite:
                        var_attrs_new[vv] = meta
                if var_attrs_new:
                    dsall_vars = create_ghcn_xarray(
                        infolder=infolder, varmeta=var_attrs_new, statmeta=station_df
                    )
                    if dsall_vars is None:
                        continue
                    dsall_vars = dsall_vars.sel(
                        time=slice(str(start_date.year), str(end_date.year))
                    )
                    for kk, vv in var_attrs_new.items():
                        cf_var = var_attrs[kk]["_cf_variable_name"]
                        outzarr = working_folder.joinpath(
                            "zarr", cf_var, f"{project}_{infolder.name}.zarr"
                        )
                        if kk not in dsall_vars.data_vars:
                            continue
                        dsout = dsall_vars.drop_vars(
                            [v for v in dsall_vars.data_vars if not v.startswith(kk)]
                        )
                        allnull_stat = dsout[kk].isnull().sum(dim="time") == len(
                            dsout.time
                        )
                        dsout = dsout.sel(station=~allnull_stat)
                        dsout = make_monotonous_time(dsout, freq=prj_dict[project])

                        ds_corr = dataset_conversion(
                            dsout,
                            project=project,
                            add_version_hashes=False,
                            overwrite=overwrite,
                        )
                        ds_corr = ds_corr.rename({f"{kk}_q_flag": f"{cf_var}_q_flag"})
                        ds_corr[f"{cf_var}_q_flag"].attrs[
                            "long_name"
                        ] = f"Quality flag for {cf_var}"
                        ds_corr[f"{cf_var}_q_flag"].attrs[
                            "description"
                        ] = f"See the readme file for information of quality flag (QFLAG1) codes : {readme_url}"
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
            if delete_raw:
                for infolder in infolders:
                    shutil.rmtree(infolder)
            treated.extend(infolders)


def make_monotonous_time(ds: xr.Dataset = None, freq: str = None):
    if freq == "daily":
        time1 = pd.date_range(ds.time[0].values, ds.time[-1].values, freq="D")
    elif freq == "hourly":
        time1 = pd.date_range(ds.time[0].values, ds.time[-1].values, freq="H")
    dsnew = xr.Dataset(coords=dict(time=time1, station=ds.station))
    for vv in ds.data_vars:
        dsnew[vv] = ds[vv]
    return dsnew


def write_zarr(
    ds: xr.Dataset = None,
    outzarr: Path = None,
    overwrite: bool = False,
    chunks: dict = None,
):
    if not outzarr.exists() or overwrite:
        with ProgressBar():
            ds.chunk(chunks).to_zarr(outzarr.with_suffix(".tmp.zarr"), mode="w")
        shutil.move(outzarr.with_suffix(".tmp.zarr"), outzarr)

        # jobs = []
        # for vv, meta in var_attrs.items():
        #     meta_var = {v[0]:v[1] for v in NoaaGhcnMetadata[prj_dict[project]].data[vv]}
        #     meta_var1 = {}
        #     for new_key, old_key in {'long_name':'name', 'units':'unit','ghcn_name':'name_original'}.items():
        #         meta_var1[new_key] = meta_var[old_key]
        #     del meta_var
        #     meta_var = meta_var1
        #     del meta_var1

        #     #jobs.append((request, ss, Path("/exec/logan/test/"), ii))

        #     create_ghcn_xarray(request=request, station_ids=ss, varmeta = meta_var, outfolder=Path("/exec/logan/test/"), izarr=ii)
