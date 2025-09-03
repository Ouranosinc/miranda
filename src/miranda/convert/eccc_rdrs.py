"""Environment and Climate Change Canada RDRS conversion tools."""

from __future__ import annotations
import logging
import os
import pathlib
from pathlib import Path
from typing import Any

import h5py
import xarray as xr
from numpy import unique

from miranda.io import fetch_chunk_config, write_dataset_dict
from miranda.treatments.utils import load_json_data_mappings
from miranda.units import check_time_frequency

from ._aggregation import aggregate
from ._data_definitions import gather_eccc_rdrs, gather_raw_rdrs_by_years
from .corrections import dataset_conversion


logger = logging.getLogger("miranda.convert.eccc_rdrs")


CONFIG_FOLDER = pathlib.Path(__file__).parent / "data"
CONFIG_FILES = {
    "rdrs-v21": "eccc_rdrs_cf_attrs.json",
    "casr-v31": "eccc_casr_cf_attrs.json",
    "ORRC-v10": "ouranos_orrc_cf_attrs.json",
    "ORRC-v11": "ouranos_orrc_cf_attrs.json",
}
for k, v in CONFIG_FILES.items():
    CONFIG_FILES[k] = CONFIG_FOLDER / v


# __all__ = ["convert_rdrs", "rdrs_to_daily"]


# FIXME: Can we use `name_output_file` instead? We already have a better version of this function.


def _get_drop_vars(file: str | os.PathLike[str], *, keep_vars: list[str] | set[str]):
    """
    Determine dropped variables.

    Parameters
    ----------
    file : str or os.PathLike
        The file to check.
    keep_vars : list or set of str
        The variables to keep.

    Returns
    -------
    list
        The dropped variables.
    """
    drop_vars = list(xr.open_dataset(file).data_vars)
    return list(set(drop_vars) - set(keep_vars))


def convert_rdrs(
    project: str,
    input_folder: str | os.PathLike[str],
    output_folder: str | os.PathLike[str],
    output_format: str = "zarr",
    working_folder: str | os.PathLike[str] | None = None,
    overwrite: bool = False,
    year_start: int | None = None,
    year_end: int | None = None,
    cfvariable_list: list | None = None,
    **dask_kwargs: dict[str, Any],
) -> None:
    r"""
    Convert RDRS dataset.

    Parameters
    ----------
    project : str
        The project name.
    input_folder : str or os.PathLike
        The input folder.
    output_folder : str or os.PathLike
        The output folder.
    output_format : {"netcdf", "zarr"}
        The output format.
    working_folder : str or os.PathLike, optional
        The working folder.
    overwrite : bool
        Whether to overwrite existing files. Default: False.
    year_start : int, optional
        The start year.
        If not provided, the minimum year in the dataset will be used.
    year_end : int, optional
        The end year.
        If not provided, the maximum year in the dataset will be used.
    cfvariable_list : list, optional
        The CF variable list.
    **dask_kwargs : dict
        Additional keyword arguments passed to the Dask scheduler.
    """
    # TODO: This setup configuration is near-universally portable. Should we consider applying it to all conversions?
    var_attrs = load_json_data_mappings(project, CONFIG_FILES)["variables"]
    prefix = load_json_data_mappings(project, CONFIG_FILES)["Header"]["_prefix"][project]
    if cfvariable_list:
        var_attrs = {
            v: var_attrs[v] for v in var_attrs if "_cf_variable_name" in var_attrs[v] and var_attrs[v]["_cf_variable_name"] in cfvariable_list
        }
    freq_dict = dict(h="hr", d="day")

    var_name_map = {v: f"{prefix}{v}" for v in var_attrs.keys()}

    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()

    # FIXME: Do we want to collect everything? Maybe return a dictionary with years and associated files?
    gathered = gather_raw_rdrs_by_years(input_folder, project)
    for year, ncfiles in gathered[project].items():
        if year_start and int(year) < year_start:
            continue
        if year_end and int(year) > year_end:
            continue
        ds_allvars = None

        if len(ncfiles) >= 28:
            for nc in ncfiles:
                eng = "h5netcdf" if h5py.is_hdf5(nc) else "netcdf4"
                try:
                    ds1 = xr.open_dataset(nc, chunks="auto", engine=eng, cache=False)
                except (OSError, RuntimeError) as e:
                    msg = f"Failed to open {nc} with engine {eng}. Error: {e}"
                    logger.error(msg)
                    raise RuntimeError(msg) from e
                if isinstance(ds1.indexes["time"], xr.coding.cftimeindex.CFTimeIndex):
                    ds1 = ds1.copy()
                    ds1["time"] = ("time", ds1.indexes["time"].to_datetimeindex())
                if ds_allvars is None:
                    out_freq = None
                    ds_allvars = ds1
                    if out_freq is None:
                        out_freq, meaning = check_time_frequency(ds1)
                        out_freq = f"{out_freq[0]}{freq_dict[out_freq[1]]}" if meaning == "hour" else freq_dict[out_freq[1]]
                    ds_allvars.attrs["frequency"] = out_freq
                else:
                    ds_allvars = xr.concat([ds_allvars, ds1], data_vars="minimal", dim="time")
            ds_allvars = ds_allvars.sel(time=f"{year}")
            # This is the heart of the conversion utility; We could apply this to multiple projects.
            for month in unique(ds_allvars.time.dt.month):
                ds_month = ds_allvars.sel(time=f"{year}-{str(month).zfill(2)}")
                for short_var in var_attrs.keys():
                    real_var_name = var_name_map[short_var]
                    drop_vars = _get_drop_vars(ncfiles[0], keep_vars=[real_var_name, "rotated_pole"])
                    ds_out = ds_month.drop_vars(drop_vars)
                    ds_out = ds_out.rename({real_var_name: short_var})
                    ds_out = ds_out.assign_coords(rotated_pole=ds_out["rotated_pole"])
                    if ds_out[short_var].attrs["units"] == "-":
                        ds_out[short_var].attrs["units"] = "1"
                    ds_corr = dataset_conversion(
                        ds_out,
                        project=project,
                        add_version_hashes=False,
                        overwrite=overwrite,
                    )
                    if "level" in ds_corr.dims:
                        ds_corr = ds_corr.squeeze()
                        ds_corr = ds_corr.drop_vars(["a", "b", "level"])
                    chunks = fetch_chunk_config(priority="time", freq=out_freq, dims=ds_corr.dims)
                    chunks["time"] = len(ds_corr.time)
                    cf_var = var_attrs[short_var]["_cf_variable_name"] if var_attrs[short_var]["_cf_variable_name"] else short_var
                    write_dataset_dict(
                        {cf_var: ds_corr},
                        output_folder=output_folder.joinpath(out_freq),
                        temp_folder=working_folder,
                        output_format=output_format,
                        overwrite=overwrite,
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
    **dask_kwargs: dict[str, Any],
) -> None:
    r"""
    Write out RDRS files to daily-timestep files.

    Parameters
    ----------
    project : str
        The project name.
    input_folder : str or os.PathLike
        The input folder.
    output_folder : str or os.PathLike
        The output folder.
    working_folder : str or os.PathLike
        The working folder.
    overwrite : bool
        Whether to overwrite existing files. Default: False.
    output_format : {"netcdf", "zarr"}
        The output format.
    year_start : int, optional
        The start year.
        If not provided, the minimum year in the dataset will be used.
    year_end : int, optional
        The end year.
        If not provided, the maximum year in the dataset will be used.
    process_variables : list of str, optional
        The variables to process.
        If not provided, all variables will be processed.
    **dask_kwargs : dict
        Additional keyword arguments passed to the Dask scheduler.
    """
    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()  # noqa
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()

    # GATHER ALL RDRS FILES
    gathered = gather_eccc_rdrs(project, input_folder, "zarr", "cf")
    files = gathered[project]  # noqa
    if process_variables:
        for vv in [f for f in files.keys() if f not in process_variables]:
            files.pop(vv)
    for zarrs in files.values():
        zarrs = sorted(zarrs)
        if not year_start:
            year_start = xr.open_zarr(zarrs[0]).time.dt.year.min().values
        if not year_end:
            year_end = xr.open_zarr(zarrs[-1]).time.dt.year.max().values
        for year in range(year_start, year_end + 1):
            infiles = [z for z in zarrs if f"_{year}" in z.name]
            if len(infiles) != 12:
                msg = f"Found {len(infiles)} input files for {year}. The year is incomplete."
                logger.warning(msg)
            out_variables = aggregate(xr.open_mfdataset(infiles, engine="zarr"), freq="day")
            dims = set(next(iter(out_variables.values())).dims)
            chunks = fetch_chunk_config(priority="time", freq="day", dims=dims)
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
