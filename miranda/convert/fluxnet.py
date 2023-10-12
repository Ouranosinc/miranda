import logging.config
import multiprocessing as mp
import os
import shutil
from pathlib import Path
from typing import List, Optional, Set, Union
from zipfile import ZipFile

import numpy as np
import pandas as pd
import xarray as xr
from dask import compute
from dask.diagnostics import ProgressBar
from dask.distributed import Client

from miranda.io.utils import delayed_write
from miranda.scripting import LOGGING_CONFIG

from ._data_corrections import dataset_conversion

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["convert_fluxnet"]


def write_job(
    ds=None,
    tmpout=None,
    out1=None,
    overwrite=None,
    output_format=None,
    target_chunks=None,
):
    if not out1.exists() or overwrite:
        tmpout.parent.mkdir(parents=True, exist_ok=True)
        job = delayed_write(
            ds,
            outfile=tmpout,
            target_chunks=target_chunks,
            output_format=output_format,
            overwrite=True,
        )
        with ProgressBar():
            compute(job)
        shutil.move(tmpout, out1)


def add_allvars(dslist: List[xr.Dataset]) -> List[xr.Dataset]:
    """Get all variables from a list of datasets.
    Parameters
    ----------
    dslist: list of xr.Dataset
    Returns
    -------
    dslist: list of xr.Dataset
    """
    list1 = [{v: d[v].attrs for v in d} for d in dslist]
    all_vars = {}
    for ll in list1:
        for key, item in ll.items():
            if key not in all_vars.keys():
                all_vars[key] = item

    for d in dslist:
        for vv in all_vars:
            if vv not in d.data_vars:
                d[vv] = xr.full_like(d["TA_F"], np.nan)
                d[vv].attrs = all_vars[vv]
    return dslist


def _combine_part(args):
    sublist, ii, tmp_folder, chunks = args
    print(f"Processing sublist {ii + 1}")
    dslist = [xr.open_dataset(f, chunks=chunks) for f in sublist]
    dslist = add_allvars(dslist)
    ds = xr.concat(dslist, dim="SITE_ID")
    tmpzarr = Path(tmp_folder) / f"tmp_{ii}.zarr"
    tmpzarr.parent.mkdir(parents=True, exist_ok=True)
    if tmpzarr.exists():
        shutil.rmtree(tmpzarr)
    with ProgressBar():
        ds.chunk(chunks).to_zarr(tmpzarr, mode="w")


def combine_station_data(
    infiles: List[Union[str, os.PathLike]],
    ds_stat: xr.Dataset,
    nsub: int = 25,
    chunks: dict = None,
    n_workers: int = 2,
    working_folder: Union[str, os.PathLike] = None,
    outpath=Union[str, os.PathLike],
    output_format: str = "zarr",
    dask_kwargs: dict = None,
):
    """Combine station data into a single file.
    Parameters
    ----------
    infiles: list of files to combine
    ds_stat: dataset with station metadata
    nsub: number of files to combine at once
    chunks: dictionary of chunks to use for each dimension
    working_folder: path to temporary working folder
    outpath: path to output file
    overwrite: whether to overwrite existing file
    output_format: output format (zarr or netcdf)
    """
    if not working_folder:
        working_folder = Path(".")
    tmpfolder = Path(working_folder).joinpath("tmpout")
    if tmpfolder.exists():
        shutil.rmtree(tmpfolder)
    tmpfolder.mkdir(parents=True, exist_ok=True)
    sublists = [infiles[x : x + nsub] for x in range(0, len(infiles), nsub)]
    for ii, sublist in enumerate(sublists):
        _combine_part((sublist, ii, tmpfolder, chunks))
    dslist = [xr.open_zarr(f, chunks=chunks) for f in tmpfolder.glob("*.zarr")]
    dslist = add_allvars(dslist)
    ds = xr.concat(dslist, dim="SITE_ID")
    # ds = xr.open_mfdataset(tmpfolder.glob("*.zarr"), engine='zarr', chunks=chunks)
    for vv in ds.data_vars:
        if "coords" in ds[vv].attrs:
            ds[vv].attrs.pop("coords")
            ds = ds.sortby(["SITE_ID", "time"])
    ds = ds_stat.sel(SITE_ID=ds.SITE_ID).merge(ds)
    ds = ds.assign_coords(
        time=pd.date_range(
            ds.time[0].values, ds.time[-1].values, freq=xr.infer_freq(dslist[0].time)
        )
    )

    return ds

    # for k,i in chunks.items():
    #     if i == -1:
    #         chunks[k] = len(ds[k])
    # tmppath = outpath.parent.joinpath(f"{outpath.stem}.tmp{outpath.suffix}")
    # tmppath.parent.mkdir(parents=True, exist_ok=True)
    # for vv in ds.coords:
    #     ds[vv] = ds[vv].load()
    # #write_job(ds=ds, tmpout=tmppath, out1=outpath, overwrite=overwrite, output_format=output_format, target_chunks=chunks)
    # job = delayed_write(ds, outfile=tmppath, target_chunks=chunks, output_format=output_format, overwrite=True)
    # with ProgressBar():
    #     compute(job)
    #
    # outpath.parent.mkdir(parents=True, exist_ok=True)
    # shutil.move(tmppath, outpath)
    # shutil.rmtree(tmpfolder)


def convert_fluxnet(
    project: str,
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    output_format: str = "zarr",
    working_folder: Optional[Union[str, os.PathLike]] = None,
    overwrite: bool = False,
    **dask_kwargs,
):
    """Convert fluxnet data to the standard format.
    Parameters
    ----------
    project : str
        Project name.
    input_folder : str or os.PathLike
        Path to the input folder.
    output_folder : str or os.PathLike
        Path to the output folder.
    output_format : str
        Output format. Default is "zarr".
    working_folder : str or os.PathLike, optional
        Path to the working folder. Default is None.
    overwrite : bool, optional
        Overwrite existing files. Default is False.
    dask_kwargs : dict
        Keyword arguments to be passed to dask.
    """
    output_format_dict = dict(zarr="zarr", netcdf="nc")
    if not working_folder:
        working_folder = Path(".")
    if not working_folder.exists():
        working_folder.mkdir(parents=True)
    data_folder = Path(__file__).parent / "data"
    var_df = pd.read_csv(
        data_folder.joinpath("fluxnet_variable_codes_SUBSET_20200504.csv")
    )
    stat_df = pd.read_csv(data_folder.joinpath("fluxnet_stations.csv"))
    ds_stat = stat_df.set_index(stat_df["SITE_ID"]).to_xarray()
    ds_stat = ds_stat.drop_vars(["MAP", "MAT", "FLUXNET-CH4"]).sel(
        SITE_ID=ds_stat["FLUXNET2015"] == "CC-BY-4.0"
    )
    ds_stat = ds_stat.rename(
        {"LOCATION_LAT": "lat", "LOCATION_LONG": "lon", "LOCATION_ELEV": "elev"}
    )
    for vv in ds_stat.data_vars:
        ds_stat = ds_stat.assign_coords({vv: ds_stat[vv]})

    for sid in ds_stat.SITE_ID.values:
        print(sid)
        zip1 = list(input_folder.glob(f"FLX_{sid}_FLUXNET2015_SUBSET_*.zip"))
        if len(zip1) != 1:
            raise ValueError(f"Found {len(zip1)} zip files for {sid}. Expected 1.")
        if zip1:
            ZipFile(zip1[0]).extractall(working_folder.joinpath(f"{sid}_csvs"))
            for freq in ["HH", "HR", "DD"]:
                tmpout = working_folder.joinpath(
                    "ind_stations",
                    freq,
                    f"{freq}_{sid}.tmp.{output_format_dict['netcdf']}",
                )
                out1 = tmpout.parent.joinpath(tmpout.name.replace(".tmp", ""))
                if not out1.exists() or overwrite:
                    freq_units = var_df[
                        var_df["Variable"].str.contains(freq.replace("HR", "HH"))
                        is True
                    ]
                    csv1 = list(
                        working_folder.joinpath(f"{sid}_csvs").glob(
                            f"FLX_{sid}_FLUXNET2015_SUBSET_{freq}_*.csv"
                        )
                    )
                    if len(csv1) > 0:
                        if len(csv1) != 1:
                            raise ValueError(
                                f"Found {len(csv1)} csv files for {sid}. Expected 1."
                            )
                        df = pd.read_csv(csv1[0])
                        time_col = (
                            "TIMESTAMP_START"
                            if "TIMESTAMP_START" in df.columns
                            else "TIMESTAMP"
                        )
                        df = df.set_index(df[time_col])
                        ds = df.to_xarray().rename({time_col: "time"})
                        if freq in ["HH", "HR"]:
                            ds = ds.assign_coords(
                                {
                                    "time": pd.to_datetime(
                                        ds.time.values, format="%Y%m%d%H%M"
                                    )
                                }
                            )
                            time1 = ds.time.values
                            time2 = xr.DataArray(
                                pd.to_datetime(
                                    ds.TIMESTAMP_END.values, format="%Y%m%d%H%M"
                                )
                            ).values
                            time_bnds = xr.DataArray(
                                list(zip(time1, time2)), dims=["time", "bnds"]
                            )

                            # ds['time_bnds'] = time_bnds
                            ds = ds.assign_coords({"time_bnds": time_bnds}).drop_vars(
                                "TIMESTAMP_END"
                            )
                            target_chunks = {"SITE_ID": 1, "time": 24 * 30}
                        elif freq == "DD":
                            target_chunks = {"SITE_ID": 1, "time": (365 * 4) + 1}
                            ds = ds.assign_coords(
                                {
                                    "time": pd.to_datetime(
                                        ds.time.values, format="%Y%m%d"
                                    )
                                }
                            )

                        for vv in ds.data_vars:
                            print(vv)
                            if vv in ["SW_OUT", "LW_OUT"]:
                                vv
                            ds[vv] = ds[vv].where(
                                ds[vv] != -9999
                            )  # convert null values to nan
                            attrs = dict()
                            if any(
                                [
                                    a in vv
                                    for a in [
                                        "TS_F_MDS_",
                                        "SWC_F_MDS_",
                                        "NEE_VUT_",
                                        "RECO_NT_VUT_",
                                        "RECO_DT_VUT_",
                                        "GPP_NT_VUT_",
                                        "GPP_DT_VUT_",
                                    ]
                                ]
                            ) and not any(
                                [
                                    a in vv
                                    for a in [
                                        "NEE_VUT_REF",
                                        "RECO_NT_VUT_REF",
                                        "RECO_DT_VUT_REF",
                                        "GPP_NT_VUT_REF",
                                        "GPP_DT_VUT_REF",
                                    ]
                                ]
                            ):
                                var1 = var_df[
                                    var_df["Variable"].str.startswith(
                                        vv.split("_QC")[0][:-2]
                                    )
                                    is True
                                ]
                                cond = var1["Variable"].str.endswith("_QC")
                                if vv.endswith("_QC"):
                                    var1 = var1[cond]
                                else:
                                    var1 = var1[~cond]
                                var1 = var1[~var1["Variable"].str.contains("_REF")]
                            else:
                                var1 = var_df[var_df["Variable"] == vv]
                            units_desc = freq_units[
                                freq_units.index.values > var1.index.values
                            ].iloc[0]
                            attrs[
                                "long_name"
                            ] = f"{var1['Description'].values[0]}".split(",")[0]
                            attrs[
                                "description"
                            ] = f"{var1['Description'].values[0]}; {units_desc['Description']}"
                            if units_desc.isnull()["Units"]:
                                attrs["units"] = ""
                            else:
                                attrs["units"] = (
                                    units_desc["Units"]
                                    .replace("deg C", "degC")
                                    .replace("nondimensional", "1")
                                )
                            ds[vv].attrs = attrs
                        ds = ds.assign_coords(
                            SITE_ID=ds_stat.sel(SITE_ID=sid).SITE_ID.values
                        ).expand_dims("SITE_ID")
                        write_job(
                            ds=ds,
                            tmpout=tmpout,
                            out1=out1,
                            output_format="netcdf",
                            overwrite=overwrite,
                        )
            shutil.rmtree(working_folder.joinpath(f"{sid}_csvs"))

    for freq, freq_out in {"DD": "day", "HH": "30min", "HR": "1hr"}.items():
        infiles = sorted(
            list(
                working_folder.joinpath(
                    "ind_stations",
                    freq,
                ).glob(f"*.{output_format_dict['netcdf']}")
            )
        )
        if len(infiles) > 0:
            if freq in ["HH", "HR"]:
                target_chunks = {"SITE_ID": 1, "time": -1}
            elif freq == "DD":
                target_chunks = {"SITE_ID": 1, "time": -1}  # (365 * 4) + 1}
            outpath = output_folder.joinpath(
                f"{freq_out}_FLUXNET2015_SUBSET_allstations.{output_format_dict[output_format]}"
            )
            if not outpath.exists() or overwrite:
                dsout = combine_station_data(
                    infiles=infiles,
                    ds_stat=ds_stat,
                    nsub=10,
                    n_workers=2,
                    chunks=target_chunks,
                    working_folder=working_folder,
                    outpath=outpath,
                    output_format=output_format,
                )
                ds_corr = dataset_conversion(
                    dsout,
                    project=project,
                    add_version_hashes=False,
                    overwrite=overwrite,
                )
                ds_corr.attrs["frequency"] = freq_out
                for k, i in target_chunks.items():
                    if i == -1:
                        target_chunks[k] = len(ds_corr[k])
                tmppath = outpath.parent.joinpath(f"{outpath.stem}.tmp{outpath.suffix}")
                tmppath.parent.mkdir(parents=True, exist_ok=True)
                job = delayed_write(
                    ds_corr,
                    outfile=tmppath,
                    target_chunks=target_chunks,
                    output_format=output_format,
                    overwrite=True,
                )
                with ProgressBar():
                    compute(job)

                outpath.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(tmppath, outpath)

    shutil.rmtree(working_folder)
