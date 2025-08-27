"""Adjusted and Homogenized Canadian Clime Data module."""

from __future__ import annotations
import calendar
import logging
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr
from dask.diagnostics import ProgressBar

from miranda.convert.utils import _add_coords_to_dataset
from miranda.eccc._utils import cf_ahccd_metadata


logger = logging.getLogger("miranda")

__all__ = [
    "convert_ahccd",
    "convert_ahccd_fwf_files",
    "create_canhomt_xarray",
    "download_canhomt",
]


def download_canhomt(
    project: str,
    working_folder: str | os.PathLike[str] | None = None,
    update_raw: bool = False,
    timeout: int | None = None,
    retry: int = 5,
    n_workers: int | None = None,  # FIXME: Not implemented yet
) -> None:
    """
    Download CanHomT data.

    Parameters
    ----------
    project : {"canhomt_dly"}
        Project name.
    working_folder : str or os.PathLink[str], optional
        Temporary files folder.
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
    if update_raw and working_folder.joinpath("raw").exists():
        shutil.rmtree(working_folder.joinpath("raw"))
    working_folder.mkdir(parents=True, exist_ok=True)

    out_folder = working_folder.joinpath("raw")
    out_folder.mkdir(parents=True, exist_ok=True)

    if project in ["canhomt_dly"]:
        url = "https://crd-data-donnees-rdc.ec.gc.ca/CDAS/products/CanHomTV4/CanHomT_dlyV4.tar.gz"
        out_file = out_folder / "CanHomT_dlyV4.tar.gz"
    else:
        msg = f"unknown project type : {project}. Supported projects are ['canhomt_dly']"
        raise ValueError(msg)

    for _ in range(retry):
        if out_file.exists() and not update_raw:
            break
        try:
            msg = f"Downloading {url}"
            logging.info(msg)
            with requests.get(url, timeout=(5, timeout)) as r:
                r.raise_for_status()
                with Path(out_file.with_suffix(f".tmp{out_file.suffix}")).open("wb") as f:
                    f.write(r.content)
            shutil.move(out_file.with_suffix(f".tmp{out_file.suffix}"), out_file)
            break
        except OSError:
            raise
    shutil.unpack_archive(out_file, out_folder, "gztar")


def create_canhomt_xarray(in_files: list, variable_meta: dict, station_meta: pd.DataFrame, project: str) -> xr.Dataset | None:
    """
    Create a xarray dataset from CanHomT raw .csv files.

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
    for cc in ["beginyear", "endyear"]:
        station_meta[cc] = pd.to_datetime(station_meta[cc], format="%Y-%m")
    for station_id in station_meta.station_id:
        msg = f"Reading {station_id}"
        logging.info(msg)
        statfile = [f for f in in_files if f.name == f"{station_id}.csv"]

        if statfile == []:
            msg = f"Station {station_id} not found in input files"
            logging.warning(msg)
            continue
        df = pd.read_csv(statfile[0])

        df.columns = df.columns.str.lower()
        df["station"] = station_id
        varlist = [k for k in variable_meta.keys() if k in df.columns]
        if varlist:
            df["time"] = pd.to_datetime(df["time"], format="%Y-%m-%d")
            df = df.set_index(["station", "time"])

            ds = df.to_xarray()
            drop_vars = [v for v in ds.data_vars if v not in varlist and "flag" not in v]
            ds = ds.drop_vars(drop_vars)
            df_stat = station_meta[station_meta.station_id == station_id]
            if len(df_stat) != 1:
                raise ValueError(f"expected a single station metadata for {station_id.stem}")
            ds = _add_coords_to_dataset(ds, df_stat)
            ds = ds.rename({f"{v}_flag": f"{v}_q_flag" for v in ds.data_vars if "flag" not in v})
            data.append(ds)

    if len(data) == 0:
        return None
    return xr.concat(data, dim="station")


def convert_ahccd(
    data_source: str | Path,
    output_dir: str | Path,
    variable: str,
    generation: int | None = None,
) -> None:
    """
    Convert Adjusted and Homogenized Canadian Climate Dataset files.

    Parameters
    ----------
    data_source: str or Path
    output_dir: str or Path
    variable: str
    generation: int, optional

    Returns
    -------
    None
    """
    output_dir = Path(output_dir).resolve().joinpath(variable)
    output_dir.mkdir(parents=True, exist_ok=True)

    code = dict(tasmax="dx", tasmin="dn", tas="dm", pr="dt", prsn="ds", prlp="dr").get(variable)
    var, col_names, col_spaces, header_row, global_attrs = cf_ahccd_metadata(code, generation)
    gen = {2: "Second", 3: "Third"}.get(generation)
    if generation == 3 and code in {"dx", "dn", "dm"}:
        meta = "ahccd_gen3_temperature.csv"
    elif generation == 2 and code in {"dt", "ds", "dr"}:
        meta = "ahccd_gen2_precipitation.csv"

    else:
        raise NotImplementedError(f"Code '{code} for generation {gen}.")
    metadata_source = Path(__file__).resolve().parent.joinpath("data").joinpath(meta)

    if "tas" in variable:
        metadata = pd.read_csv(metadata_source, header=2)
        metadata.columns = col_names.keys()
        cols_specs = col_spaces

    elif "pr" in variable:
        metadata = pd.read_csv(metadata_source, header=3)
        metadata.columns = col_names.keys()
        cols_specs = col_spaces
        for index, row in metadata.iterrows():
            if isinstance(row["stnid"], str):
                metadata.loc[index, "stnid"] = metadata.loc[index, "stnid"].replace(" ", "")
    else:
        raise KeyError(f"{variable} does not include 'pr' or 'tas'.")

    # Convert station .txt files to netcdf
    for ff in Path(data_source).glob("*d*.txt"):
        outfile = output_dir.joinpath(ff.name.replace(".txt", ".nc"))
        if not outfile.exists():
            logger.info(ff.name)

            stid = ff.name.replace(code, "").split(".txt")[0]
            try:
                metadata_st = metadata[metadata["stnid"] == int(stid)]
            except ValueError:
                metadata_st = metadata[metadata["stnid"] == stid]

            if len(metadata_st) == 1:
                ds_out = convert_ahccd_fwf_files(ff, metadata_st, variable, generation, cols_specs, var)
                ds_out.attrs = global_attrs

                ds_out.to_netcdf(outfile, engine="h5netcdf")
            else:
                msg = f"metadata info for station {ff.name} not found : skipping"

                logger.warning(msg)

    # merge individual stations to single .nc file
    # variable
    ncfiles = list(output_dir.glob("*.nc"))
    outfile = output_dir.parent.joinpath("merged_stations", f"ahccd_gen{generation}_{variable}.nc")

    if not outfile.exists():
        logger.info("merging stations :", variable)
        with ProgressBar():
            ds_ahccd = xr.open_mfdataset(ncfiles, concat_dim="station", combine="nested").load()

            for coord in ds_ahccd.coords:
                # xarray object datatypes mix string and int (e.g. stnid) convert to string for merged nc files
                # Do not apply to datetime object
                if coord != "time" and ds_ahccd[coord].dtype == "O":
                    ds_ahccd[coord] = ds_ahccd[coord].astype(str)

            for v in ds_ahccd.data_vars:
                # xarray object datatypes mix string and int (e.g. stnid) convert to string for merged nc files
                # Do not apply to flag timeseries
                if ds_ahccd[v].dtype == "O" and "flag" not in v:
                    logger.info(v)
                    ds_ahccd[v] = ds_ahccd[v].astype(str)

            ds_ahccd[f"{variable}_flag"].attrs["long_name"] = f"{ds_ahccd[f'{variable}'].attrs['long_name']} flag"
            ds_ahccd.lon.attrs["units"] = "degrees_east"
            ds_ahccd.lon.attrs["long_name"] = "longitude"
            ds_ahccd.lat.attrs["units"] = "degrees_north"
            ds_ahccd.lat.attrs["long_name"] = "latitude"

            for clean_name, orig_name in col_names.items():
                if clean_name in ["lat", "long"]:
                    continue
                ds_ahccd[clean_name].attrs["long_name"] = orig_name

            outfile.parent.mkdir(parents=True, exist_ok=True)
            ds_ahccd.to_netcdf(outfile, engine="h5netcdf", format="NETCDF4_CLASSIC", mode="w")

            del ds_ahccd
    for nc in outfile.parent.glob("*.nc"):
        logger.info(nc)
        ds = xr.open_dataset(nc)
        logger.info(ds)


def convert_ahccd_fwf_files(
    ff: Path | str,
    metadata: pd.DataFrame,
    variable: str,
    generation: int | None = None,
    cols_specs: list[tuple[int, int]] | None = None,
    attrs: dict | None = None,
) -> xr.Dataset:
    """
    Convert AHCCD fixed-width files.

    Parameters
    ----------
    ff: str or Path
    metadata: pandas.DataFrame
    variable: str
    generation
    cols_specs
    attrs

    Returns
    -------
    xarray.Dataset
    """
    code = dict(tasmax="dx", tasmin="dn", tas="dm", pr="dt", prsn="ds", prlp="dr").get(variable)

    if attrs is None:
        attrs, _, _, _, _ = cf_ahccd_metadata(code, generation)
    if cols_specs is None:
        _, _, cols_specs, _, _ = cf_ahccd_metadata(code, generation)
    _, _, _, nhead, _ = cf_ahccd_metadata(code, generation)

    df = pd.read_fwf(ff, header=nhead, colspecs=cols_specs)
    if "pr" in variable:
        cols = list(df.columns[0:3])
        cols = cols[0::2]
        cols.extend(list(df.columns[4::2]))
        flags = list(df.columns[5::2])
        dfflags = df[flags]
    else:
        cols = [c for c in df.columns if "Unnamed" not in c]
        flags = [c for c in df.columns if "Unnamed" in c]
        dfflags = df[flags[2:]]

    df = df[cols]
    df.replace(attrs["NaN_value"], np.NaN, inplace=True)

    for i, j in enumerate(["Year", "Month"]):
        df = df.rename(columns={df.columns[i]: j})
    start_date = f"{df['Year'][0]}-{str(df['Month'][0]).zfill(2)}-01"

    _, ndays = calendar.monthrange(df["Year"].iloc[-1], df["Month"].iloc[-1])
    end_date = f"{df['Year'].iloc[-1]}-{str(df['Month'].iloc[-1]).zfill(2)}-{str(ndays).zfill(2)}"
    time1 = pd.date_range(start=start_date, end=end_date)

    index = pd.MultiIndex.from_arrays([df["Year"], df["Month"]])
    df.index = index
    dfflags.index = index
    cols = [c for c in df.columns if "Year" not in c and "Month" not in c]
    df = df[cols]
    df.columns = np.arange(1, 32)
    dfflags.columns = np.arange(1, 32)
    ds = df.stack().to_frame()
    ds = ds.rename(columns={0: variable})
    ds_flag = dfflags.stack().to_frame()
    ds_flag = ds_flag.rename(columns={0: "flag"})
    ds.index.names = ["Year", "Month", "Day"]
    ds_flag.index.names = ["Year", "Month", "Day"]
    ds[f"{variable}_flag"] = ds_flag["flag"]
    del ds_flag

    # find invalid dates
    for y in time1.year.unique():
        for m in ds[ds.index.get_level_values("Year") == y].index.get_level_values("Month").unique():
            _, exp_ndays = calendar.monthrange(y, m)
            ndays = ((ds.index.get_level_values("Year") == y) & (ds.index.get_level_values("Month") == m)).sum()
            if ndays > np.int(exp_ndays):
                print(f"year {y}, month {m}, ndays={ndays}, exp_ndays={exp_ndays}")
                raise RuntimeError("Unknown days present.")

    time_ds = pd.DataFrame(
        {
            "year": ds.index.get_level_values("Year"),
            "month": ds.index.get_level_values("Month"),
            "day": ds.index.get_level_values("Day"),
        }
    )

    ds.index = pd.to_datetime(time_ds)

    ds = ds.to_xarray().rename({"index": "time"})

    ds_out = xr.Dataset(coords={"time": time1})
    for v in ds.data_vars:
        ds_out[v] = ds[v]

    ds_out[variable].attrs = attrs
    # ds_out
    metadata = metadata.to_xarray().rename({"index": "station"}).drop_vars("station")
    metadata = metadata.assign_coords(
        {
            "stnid": metadata["stnid"].astype(str),
            "station_name": metadata["station_name"],
        }
    )
    # ds_out = ds_out.assign_coords({'lon': metadata['long'], 'lat': metadata['lat'], 'elevation': metadata['elev']})
    #
    ds_out = ds_out.assign_coords(station=metadata.stnid)
    metadata = metadata.drop_vars(["stnid", "station_name"])

    ds_out["lon"] = metadata["long"]
    ds_out["lon"].attrs["units"] = "degrees_east"
    ds_out["lat"] = metadata["lat"]
    ds_out["lat"].attrs["units"] = "degrees_north"
    ds_out["elev"] = metadata["elev"]
    ds_out["elev"].attrs["units"] = "m"

    metadata = metadata.drop_vars(["long", "lat", "elev"])
    for vv in metadata.data_vars:
        if metadata[vv].dtype == "O" and (variable not in vv):
            ds_out[vv] = metadata[vv].astype(str)
        else:
            ds_out[vv] = metadata[vv]
    return ds_out
