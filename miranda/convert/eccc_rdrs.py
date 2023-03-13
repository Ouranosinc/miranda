import datetime
import logging.config
import os
import shutil
from pathlib import Path
from typing import List, Set, Union

import xarray as xr
from dask import compute
from dask.distributed import Client

from miranda.convert import file_conversion, load_json_data_mappings
from miranda.convert.utils import daily_aggregation, delayed_write
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["convert_rdrs", "concat_zarr", "rdrs_to_daily"]


def _rename_output_file(
    fname: str, old_variable: str, new_variable: str, freq: str, date: Union[str, int]
):
    parts = fname.replace(old_variable, new_variable).split("_")
    parts[-1] = f"{date}"
    parts[1] = freq
    return f"{'_'.join(parts)}.zarr"


def _get_drop_vars(
    file: Union[str, os.PathLike], *, keep_vars: Union[List[str], Set[str]]
):
    drop_vars = list(xr.open_dataset(file).data_vars)
    return list(set(drop_vars) - set(keep_vars))


def convert_rdrs(
    project: str,
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    output_format: str = "zarr",
    working_folder: Union[str, os.PathLike] = None,
    overwrite: bool = False,
    **dask_kwargs,
):
    var_attrs = load_json_data_mappings(project=project)["variable_entry"]
    freq_dict = dict(h="hr", d="day")

    if isinstance(input_folder, str):
        input_folder = Path(input_folder)
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)
    if isinstance(working_folder, str):
        working_folder = Path(working_folder)

    if working_folder:
        if working_folder.exists():
            if overwrite:
                shutil.rmtree(working_folder)
            else:
                raise FileExistsError(
                    "`working_folder` is not empty. Use `overwrite=True`."
                )
        working_folder.mkdir(parents=True)

    year = datetime.date.today().year
    for year in range(1950, year):
        if not list(input_folder.glob(f"{year}*.nc")):
            logging.warning(f"no files found for year {year}. Continuing...")
            continue
        # Time stamps starts at noon and flow into subsequent months!!
        # Need full year plus previous december in order to easily produce complete hourly frequency monthly files
        ncfiles = sorted(list(input_folder.glob(f"{year - 1}12*.nc")))
        if ncfiles:
            ncfiles = [ncfiles[-1]]
        ncfiles.extend(sorted(list(input_folder.glob(f"{year}*.nc"))))
        ds_allvars = None
        if ncfiles:
            for nc in ncfiles:
                ds1 = xr.open_dataset(nc, chunks="auto")
                if ds_allvars is None:
                    ds_allvars = ds1
                    outfreq, meaning = get_time_frequency(ds1)
                    outfreq = (
                        f"{outfreq[0]}{freq_dict[outfreq[1]]}"
                        if meaning == "hour"
                        else freq_dict[outfreq[1]]
                    )

                else:
                    ds_allvars = xr.concat([ds_allvars, ds1], dim="time")
            ds_allvars.attrs["Remarks"] = ds_allvars.attrs["Remarks"].replace(
                "Variable names", "Original variable names"
            )
            for month in range(1, 13):
                ds_month = ds_allvars.sel(time=f"{year}-{str(month).zfill(2)}")
                for var_attr in var_attrs.keys():
                    drop_vars = _get_drop_vars(
                        ncfiles[0], keep_vars=[var_attr, "rotated_pole"]
                    )

                    # outfile_root = ncfiles[0].stem.split()
                    ds_out = ds_month.drop_vars(drop_vars)

                    overwrite_flag = False
                    write_folder = (
                        working_folder if working_folder is not None else output_folder
                    )

                    job = file_conversion(
                        input=ds_out,
                        project=project,
                        output_path=write_folder.joinpath(
                            var_attrs[var_attr]["_cf_variable_name"]
                        ),
                        output_format=output_format,
                        add_version_hashes=False,
                        chunks=dict(time=len(ds_out.time), rlon=50, rlat=50),
                        compute=False,
                        correct_variable_names=True,
                        overwrite=True,
                    )

                    if working_folder:
                        for subdir in [
                            f for f in working_folder.glob("*") if f.is_dir()
                        ]:
                            output_folder.joinpath(outfreq, subdir.name).mkdir(
                                exist_ok=True
                            )
                            for zarr in subdir.glob("*.zarr"):
                                new_zarr = zarr.name.split("_")
                                ## TODO manage output file name
                                new_zarr = f"{new_zarr[0]}_{outfreq}_eccc_{new_zarr[1]}_NAM_{new_zarr[2]}"

                                if (
                                    overwrite_flag
                                    or not output_folder.joinpath(
                                        outfreq, subdir.name, new_zarr
                                    ).exists()
                                ):
                                    with Client(**dask_kwargs):
                                        job.compute()
                                    shutil.copytree(
                                        zarr,
                                        output_folder.joinpath(
                                            outfreq, subdir.name, new_zarr.name
                                        ),
                                        dirs_exist_ok=True,
                                    )
                                else:
                                    logging.warning(
                                        output_folder.joinpath(
                                            outfreq, subdir.name, new_zarr
                                        ),
                                        " exists. Continuing...",
                                    )
                            shutil.rmtree(subdir)


def rdrs_to_daily(
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    working_folder: Union[str, os.PathLike] = None,
    overwrite: bool = False,
    **dask_kwargs,
):
    freq = "day"
    if isinstance(input_folder, str):
        input_folder = Path(input_folder)
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

    for var_folder in [v for v in input_folder.glob("*") if v.is_dir()]:
        zarrs = sorted(list(var_folder.glob("*.zarr")))
        year1 = xr.open_zarr(zarrs[0]).time.dt.year.min().values
        year2 = xr.open_zarr(zarrs[-1]).time.dt.year.max().values
        for year in range(year1, year2 + 1):
            infiles = sorted(list(var_folder.glob(f"*_{year}*.zarr")))

            if len(infiles) > 12 or len(infiles) == 0:
                raise ValueError(
                    f"found {len(infiles)} input files ... expected 12 or less"
                )
            outvars = daily_aggregation(xr.open_zarr(infiles[0]), keys_only=True)
            for vv in outvars:
                ds = None
                outzarr = _rename_output_file(
                    infiles[0].stem, var_folder.name, vv, freq, year
                )
                outzarr = output_folder.joinpath(vv, outzarr)
                if not outzarr.exists() or overwrite:
                    if ds is None:
                        ds = xr.open_mfdataset(infiles, engine="zarr")
                        ds = daily_aggregation(ds=ds)
                    outzarr.parent.mkdir(parents=True, exist_ok=True)
                    if outzarr.exists():
                        shutil.rmtree(outzarr)

                    if working_folder:
                        tmp_zarr = working_folder.joinpath(outzarr.name)
                        if tmp_zarr.exists():
                            shutil.rmtree(tmp_zarr)
                        with Client(**dask_kwargs):
                            ds[vv].chunk(dict(time=-1, rlon=50, rlat=50)).to_zarr(
                                tmp_zarr
                            )
                        shutil.copytree(tmp_zarr, outzarr, dirs_exist_ok=True)
                        shutil.rmtree(tmp_zarr)
                    else:
                        with Client(**dask_kwargs):
                            ds[var_folder.name].chunk(
                                dict(time=-1, rlon=50, rlat=50)
                            ).to_zarr(outzarr)
                else:
                    print(outzarr.as_posix(), "exists. Continuing..")


def concat_zarr(
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    overwrite: bool = False,
    **dask_kwargs,
):
    if isinstance(input_folder, str):
        input_folder = Path(input_folder)
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

    list_zarr = sorted(list(input_folder.glob("*.zarr")))
    outzarr = "_".join(list_zarr[0].stem.split("_")[0:-1])
    st_yr = list_zarr[0].stem.split("_")[-1].split("-")[0][0:4]
    end_yr = list_zarr[-1].stem.split("_")[-1].split("-")[0][0:4]
    outzarr = f"{outzarr}_{st_yr}_{end_yr}.zarr"
    outzarr = output_folder.joinpath(outzarr)

    if not outzarr.exists() or overwrite:
        if "day" in input_folder.as_posix():
            chunks = dict(time=(365 * 4) + 1, rlon=50, rlat=50)
        else:
            chunks = dict(time=(24 * 30 * 2), rlon=50, rlat=50)

        # maketemp files 1 zarr per 4 years
        years = [y for y in range(int(st_yr), int(end_yr) + 1)]
        years = [years[x : x + 4] for x in range(0, len(years), 4)]
        for year in years:
            list_zarr1 = sorted(
                [
                    zarrfile
                    for zarrfile in list_zarr
                    if int(zarrfile.stem.split("_")[-1].split("-")[0][0:4]) in year
                ]
            )
            assert len(list_zarr1) / len(year) == 12
            ds = xr.open_mfdataset(list_zarr1, parallel=True, engine="zarr")

            with Client(**dask_kwargs):
                # if outzarr.exists():
                #     zarr_kwargs = {"append_dim": "time", "consolidated": True}
                # else:
                #     zarr_kwargs = {"consolidated": True}
                tmpzarr = outzarr.parent.joinpath(
                    "tmp",
                    f"{outzarr.stem.split(f'_{st_yr}_')[0]}_{year[0]}-{year[-1]}.zarr",
                )
                tmpzarr.parent.mkdir(exist_ok=True, parents=True)
                logging.info(f"{year} writing to {tmpzarr.as_posix()}")

                job = delayed_write(
                    ds=ds,
                    outfile=tmpzarr,
                    output_format="zarr",
                    target_chunks=chunks,
                    overwrite=overwrite,
                )  # kwargs=zarr_kwargs)
                compute(job)

        # get tmp zarrs
        list_zarr = sorted(list(tmpzarr.parent.glob("*zarr")))
        ds = xr.open_mfdataset(list_zarr, engine="zarr")
        with Client(**dask_kwargs):
            job = delayed_write(
                ds=ds,
                outfile=outzarr,
                output_format="zarr",
                target_chunks=chunks,
                overwrite=overwrite,
            )  # kwargs=zarr_kwargs)
            compute(job)

        shutil.rmtree(tmpzarr.parent)
