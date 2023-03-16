import datetime
import logging.config
import os
import shutil
from pathlib import Path
from typing import List, Optional, Set, Union

import xarray as xr
from dask.distributed import Client

from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

from ._data_corrections import dataset_conversion, load_json_data_mappings
from ._data_definitions import gather_raw_rdrs_by_years, gather_rdrs
from ..io import write_dataset, fetch_chunks
from ..io._output import write_dataset_dict

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["convert_rdrs", "rdrs_to_daily"]


# FIXME: Can we use `name_output_file` instead? We already have a better version of this function.


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
    working_folder: Optional[Union[str, os.PathLike]] = None,
    overwrite: bool = False,
    **dask_kwargs,
):
    """

    Parameters
    ----------
    project
    input_folder
    output_folder
    output_format
    working_folder
    overwrite
    dask_kwargs

    Returns
    -------

    """
    # TODO: This setup configuration is near-universally portable. Should we consider applying it to all conversions?
    var_attrs = load_json_data_mappings(project=project)["variable_entry"]
    freq_dict = dict(h="hr", d="day")

    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()



    # FIXME: Do we want to collect everything? Maybe return a dictionary with years and associated files?

    gathered = gather_raw_rdrs_by_years(input_folder)
    for year, ncfiles in gathered[project].items():
        ds_allvars = None
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

        # This is the heart of the conversion utility; We could apply this to multiple projects.
        for month in range(1, 13):
            ds_month = ds_allvars.sel(time=f"{year}-{str(month).zfill(2)}")
            for var_attr in var_attrs.keys():
                drop_vars = _get_drop_vars(
                    ncfiles[0], keep_vars=[var_attr, "rotated_pole"]
                )
                ds_out = ds_month.drop_vars(drop_vars)


                ds_corr = dataset_conversion(
                    ds_out,
                    project=project,
                    add_version_hashes=False,
                    overwrite=overwrite,
                )
                chunks = fetch_chunks(project=project, freq="1hr")
                chunks["time"] = len(ds_corr.time)
                write_dataset_dict(
                    {var_attr:ds_corr},
                    output_folder=output_folder,
                    temp_folder=working_folder,
                    output_format=output_format,
                    overwrite=False,
                    chunks=chunks,
                    **dask_kwargs
                )

                # write_dataset(
                #     ds_corr,
                #     project=project,
                #     output_path=write_folder,
                #     output_format=output_format,
                # )

                #
                #
                # # FIXME: This is generalized IO code - should be factored out.
                # if working_folder:
                #     for subdir in [f for f in working_folder.glob("*") if f.is_dir()]:
                #         output_folder.joinpath(outfreq, subdir.name).mkdir(
                #             exist_ok=True
                #         )
                #         for zarr in subdir.glob("*.zarr"):
                #             new_zarr = zarr.name.split("_")
                #             # TODO: manage output file name
                #             new_zarr = f"{new_zarr[0]}_{outfreq}_eccc_{new_zarr[1]}_NAM_{new_zarr[2]}"
                #
                #             if (
                #                 overwrite
                #                 or not output_folder.joinpath(
                #                     outfreq, subdir.name, new_zarr
                #                 ).exists()
                #             ):
                #                 # FIXME: Client is only needed for computation. Should be elsewhere.
                #                 with Client(**dask_kwargs):
                #                     job.compute()
                #                 shutil.copytree(
                #                     zarr,
                #                     output_folder.joinpath(
                #                         outfreq, subdir.name, new_zarr.name
                #                     ),
                #                     dirs_exist_ok=True,
                #                 )
                #             else:
                #                 logging.warning(
                #                     output_folder.joinpath(
                #                         outfreq, subdir.name, new_zarr
                #                     ),
                #                     " exists. Continuing...",
                #                 )
                #         shutil.rmtree(subdir)


# FIXME: This looks mostly like code to stage writing out files. Should it be moved to an IO module?
def rdrs_to_daily(
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    working_folder: Optional[Union[str, os.PathLike]] = None,
    overwrite: bool = False,
    **dask_kwargs,
):
    """

    Parameters
    ----------
    input_folder
    output_folder
    working_folder
    overwrite
    dask_kwargs

    Returns
    -------

    """
    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()  # noqa
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()

    if working_folder:
        if working_folder.exists():
            if overwrite:
                shutil.rmtree(working_folder)
            else:
                raise FileExistsError(
                    "`temp_folder` is not empty. Use `overwrite=True`."
                )
        working_folder.mkdir(parents=True)

    # GATHER ALL RDRS FILES
    gathered = gather_rdrs(input_folder, "zarr")
    files = gathered["rdrs-v2.1"]  # noqa

    # for var_folder in [v for v in input_folder.glob("*") if v.is_dir()]:
    #     zarrs = sorted(list(var_folder.glob("*.zarr")))
    #     year1 = xr.open_zarr(zarrs[0]).time.dt.year.min().values
    #     year2 = xr.open_zarr(zarrs[-1]).time.dt.year.max().values
    #     for year in range(year1, year2 + 1):
    #         infiles = sorted(list(var_folder.glob(f"*_{year}*.zarr")))
    #         if len(infiles) > 12 or len(infiles) == 0:
    #             raise ValueError(
    #                 f"Found {len(infiles)} input files. Expected between 1 and 12."
    #             )
    #
    #         out_variables = daily_aggregation(xr.open_zarr(infiles[0]), keys_only=True)
    #         for variable in out_variables:
    #             ds = None
    #             outzarr = _rename_output_file(
    #                 infiles[0].stem, var_folder.name, variable, "day", year
    #             )
    #             outzarr = output_folder.joinpath(variable, outzarr)
    #             if not outzarr.exists() or overwrite:
    #                 if ds is None:
    #                     ds = xr.open_mfdataset(infiles, engine="zarr")
    #                     ds = daily_aggregation(ds=ds)
    #                 outzarr.parent.mkdir(parents=True, exist_ok=True)
    #                 if outzarr.exists():
    #                     shutil.rmtree(outzarr)
    #
    #                 if temp_folder:
    #                     tmp_zarr = temp_folder.joinpath(outzarr.name)
    #                     # FIXME: Client is only needed for computation. Should be elsewhere.
    #                     with Client(**dask_kwargs):
    #                         ds[variable].chunk(dict(time=-1, rlon=50, rlat=50)).to_zarr(
    #                             tmp_zarr
    #                         )
    #                     shutil.copytree(tmp_zarr, outzarr, dirs_exist_ok=True)
    #                     shutil.rmtree(tmp_zarr)
    #                 else:
    #                     with Client(**dask_kwargs):
    #                         ds[var_folder.name].chunk(
    #                             dict(time=-1, rlon=50, rlat=50)
    #                         ).to_zarr(outzarr)
    #             else:
    #                 logging.warning(outzarr.as_posix(), "exists. Continuing..")
