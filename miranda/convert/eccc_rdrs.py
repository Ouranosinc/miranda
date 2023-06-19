"""Environment and Climate Change Canada RDRS conversion tools."""
from __future__ import annotations

import logging.config
import os
from pathlib import Path

import xarray as xr
from numpy import unique

from miranda.io import fetch_chunk_config, write_dataset_dict
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

from ._aggregation import aggregate
from ._data_corrections import dataset_conversion, load_json_data_mappings
from ._data_definitions import gather_raw_rdrs_by_years, gather_rdrs

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["convert_rdrs", "rdrs_to_daily"]


# FIXME: Can we use `name_output_file` instead? We already have a better version of this function.


def _get_drop_vars(file: str | os.PathLike, *, keep_vars: list[str] | set[str]):
    drop_vars = list(xr.open_dataset(file).data_vars)
    return list(set(drop_vars) - set(keep_vars))


def convert_rdrs(
    project: str,
    input_folder: str | os.PathLike,
    output_folder: str | os.PathLike,
    output_format: str = "zarr",
    working_folder: str | os.PathLike | None = None,
    overwrite: bool = False,
    **dask_kwargs,
) -> None:
    r"""

    Parameters
    ----------
    project : str
    input_folder : str or os.PathLike
    output_folder : str or os.PathLike
    output_format : {"netcdf", "zarr"}
    working_folder : str or os.PathLike, optional
    overwrite : bool
    \*\*dask_kwargs

    Returns
    -------
    None
    """
    # TODO: This setup configuration is near-universally portable. Should we consider applying it to all conversions?
    var_attrs = load_json_data_mappings(project=project)["variables"]
    freq_dict = dict(h="hr", d="day")

    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()

    # FIXME: Do we want to collect everything? Maybe return a dictionary with years and associated files?

    out_freq = None
    gathered = gather_raw_rdrs_by_years(input_folder)
    for year, ncfiles in gathered[project].items():
        ds_allvars = None
        if len(ncfiles) >= 28:
            for nc in ncfiles:
                ds1 = xr.open_dataset(nc, chunks="auto")
                if ds_allvars is None and out_freq is None:
                    ds_allvars = ds1
                    out_freq, meaning = get_time_frequency(ds1)
                    out_freq = (
                        f"{out_freq[0]}{freq_dict[out_freq[1]]}"
                        if meaning == "hour"
                        else freq_dict[out_freq[1]]
                    )
                    ds_allvars.attrs["frequency"] = out_freq
                else:
                    ds_allvars = xr.concat(
                        [ds_allvars, ds1], data_vars="minimal", dim="time"
                    )
            ds_allvars = ds_allvars.sel(time=f"{year}")
            # This is the heart of the conversion utility; We could apply this to multiple projects.
            for month in unique(ds_allvars.time.dt.month):
                ds_month = ds_allvars.sel(time=f"{year}-{str(month).zfill(2)}")
                for var_attr in var_attrs.keys():
                    drop_vars = _get_drop_vars(
                        ncfiles[0], keep_vars=[var_attr, "rotated_pole"]
                    )
                    ds_out = ds_month.drop_vars(drop_vars)
                    ds_out = ds_out.assign_coords(rotated_pole=ds_out["rotated_pole"])
                    ds_corr = dataset_conversion(
                        ds_out,
                        project=project,
                        add_version_hashes=False,
                        overwrite=overwrite,
                    )
                    chunks = fetch_chunk_config(
                        priority="time", freq=out_freq, dims=ds_corr.dims
                    )
                    chunks["time"] = len(ds_corr.time)
                    write_dataset_dict(
                        {var_attrs[var_attr]["_cf_variable_name"]: ds_corr},
                        output_folder=output_folder.joinpath(out_freq),
                        temp_folder=working_folder,
                        output_format=output_format,
                        overwrite=False,
                        chunks=chunks,
                        **dask_kwargs,
                    )


# FIXME: This looks mostly like code to stage writing out files. Should it be moved to an IO module?
def rdrs_to_daily(
    project: str,
    input_folder: str | os.PathLike,
    output_folder: str | os.PathLike,
    working_folder: str | os.PathLike | None = None,
    overwrite: bool = False,
    output_format: str = "zarr",
    year_start: int | None = None,
    year_end: int | None = None,
    process_variables: list[str] | None = None,
    **dask_kwargs,
) -> None:
    r"""Write out RDRS files to daily-timestep files.

    Parameters
    ----------
    project : str
    input_folder : str or os.PathLike
    output_folder : str or os.PathLike
    working_folder : str or os.PathLike
    overwrite : bool
    output_format : {"netcdf", "zarr"}
    year_start : int, optional
    year_end : int, optional
    process_variables : list of str, optional
    \*\*dask_kwargs

    Returns
    -------
    None
    """
    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()  # noqa
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()

    # GATHER ALL RDRS FILES
    gathered = gather_rdrs(project, input_folder, "zarr", "cf")
    files = gathered["rdrs-v21"]  # noqa
    if process_variables:
        for vv in [f for f in files.keys() if f not in process_variables]:
            files.pop(vv)
    for vv, zarrs in files.items():
        zarrs = sorted(zarrs)
        if not year_start:
            year_start = xr.open_zarr(zarrs[0]).time.dt.year.min().values
        if not year_end:
            year_end = xr.open_zarr(zarrs[-1]).time.dt.year.max().values
        for year in range(year_start, year_end + 1):
            infiles = [z for z in zarrs if f"_{year}" in z.name]
            if len(infiles) != 12:
                raise ValueError(f"Found {len(infiles)} input files. Expected 12.")
            #
            out_variables = aggregate(
                xr.open_mfdataset(infiles, engine="zarr"), freq="day"
            )
            chunks = fetch_chunk_config(project=project, freq="day")
            chunks["time"] = len(out_variables[list(out_variables.keys())[0]].time)
            write_dataset_dict(
                out_variables,
                output_folder=output_folder,
                temp_folder=working_folder,
                output_format=output_format,
                overwrite=overwrite,
                chunks=chunks,
                **dask_kwargs,
            )
