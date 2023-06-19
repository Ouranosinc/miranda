"""IO Utilities module."""
from __future__ import annotations

import json
import logging.config
import os
from collections.abc import Sequence
from datetime import date
from pathlib import Path

import dask
import netCDF4 as nc  # noqa
import xarray as xr
import zarr

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = [
    "creation_date",
    "delayed_write",
    "get_chunks_on_disk",
    "get_global_attrs",
    "get_time_attrs",
    "name_output_file",
    "sort_variables",
]

_data_folder = Path(__file__).parent / "data"
name_configurations = json.load(open(_data_folder / "ouranos_name_config.json"))


def name_output_file(
    ds_or_dict: xr.Dataset | dict[str, str], output_format: str
) -> str:
    """Name an output file based on facets within a Dataset or a dictionary.

    Parameters
    ----------
    ds_or_dict : xr.Dataset or dict
        A miranda-converted Dataset or a dictionary containing the appropriate facets.
    output_format : {"netcdf", "zarr"}
        Output filetype to be used for generating filename suffix.

    Returns
    -------
    str

    Notes
    -----
    If using a dictionary, the following keys must be set:
    * "variable", "frequency", "institution", "time_start", "time_end".
    """
    if output_format.lower() not in {"netcdf", "zarr"}:
        raise NotImplementedError(f"Format: {output_format}.")
    else:
        suffix = dict(netcdf="nc", zarr="zarr")[output_format]

    facets = dict()
    facets["suffix"] = suffix

    if isinstance(ds_or_dict, xr.Dataset):
        if len(ds_or_dict.data_vars) == 1:
            facets["variable"] = list(ds_or_dict.data_vars.keys())[0]
        elif (
            len(ds_or_dict.data_vars) == 2
            and "rotated_pole" in ds_or_dict.data_vars.keys()
        ):
            facets["variable"] = [
                v for v in ds_or_dict.data_vars if v != "rotated_pole"
            ][0]
        else:
            raise NotImplementedError(
                f"Too many `data_vars` in Dataset: {' ,'.join(ds_or_dict.data_vars.keys())}."
            )
        for f in [
            "bias_adjust_project",
            "domain",
            "frequency",
            "institution",
            "source",
            "experiment",
            "member",
            "processing_level",
            "project",
            "type",
            "mip_era",
            "activity",
        ]:
            facets[f] = ds_or_dict.attrs.get(f)

        if facets["frequency"] in ["1hr", "day"]:
            date_format = "%Y%m%d"
        elif facets["frequency"] == "month":
            date_format = "%Y%m"
        elif facets["frequency"] == "year":
            date_format = "%Y"
        else:
            raise KeyError("`frequency` not found.")

        facets["time_start"], facets["time_end"] = (
            ds_or_dict.time.isel(time=[0, -1]).dt.strftime(date_format).values
        )
        facets["year_start"], facets["year_end"] = ds_or_dict.time.isel(
            time=[0, -1]
        ).dt.year.values
    elif isinstance(ds_or_dict, dict):
        for f in [
            "bias_adjust_project",
            "domain",
            "frequency",
            "institution",
            "processing_level",
            "project",
            "type",
            "time",
            "time_end",
            "time_start",
            "variable",
        ]:
            facets[f] = ds_or_dict.get(f)
    else:
        raise NotImplementedError("Must be a Dataset or dictionary.")

    if {"time_start", "time_end"}.issubset(facets) and "time" not in facets:
        if facets["time_start"] == facets["time_end"]:
            facets["time"] = "-".join([facets["time_start"], facets["time_end"]])
        else:
            facets["time"] = facets["time_start"]

    str_name = "{variable}_{frequency}_{institution}_{project}_{time}.{suffix}"
    # Get the string for the name
    if facets["type"] in name_configurations.keys():
        if facets["project"] in name_configurations[facets["type"]].keys():
            str_name = name_configurations[facets["type"]][facets["project"]]

    missing = []
    for k, v in facets.items():
        if (
            v is None and k in str_name
        ):  # only missing if the facets is needed in the name
            missing.append(k)
    if missing:
        raise ValueError(f"The following facets were not found: {' ,'.join(missing)}.")

    # fill in string with facets
    return str_name.format(**facets)


def delayed_write(
    ds: xr.Dataset,
    outfile: str | os.PathLike,
    output_format: str,
    overwrite: bool,
    target_chunks: dict | None = None,
) -> dask.delayed:
    """Stage a Dataset writing job using `dask.delayed` objects.

    Parameters
    ----------
    ds : xr.Dataset
    outfile : str or os.PathLike
    target_chunks : dict
    output_format : {"netcdf", "zarr"}
    overwrite : bool

    Returns
    -------
    dask.delayed.delayed
    """
    # Set correct chunks in encoding options
    kwargs = dict()
    kwargs["encoding"] = dict()
    try:
        for name, da in ds.data_vars.items():
            chunks = list()
            for dim in da.dims:
                if target_chunks:
                    if dim in target_chunks.keys():
                        chunks.append(target_chunks[str(dim)])
                else:
                    chunks.append(len(da[dim]))

            if output_format == "netcdf":
                kwargs["encoding"][name] = {
                    "chunksizes": chunks,
                    "zlib": True,
                }
                kwargs["compute"] = False
                if Path(outfile).exists() and not overwrite:
                    kwargs["mode"] = "a"
            elif output_format == "zarr":
                ds = ds.chunk(target_chunks)
                kwargs["encoding"][name] = {
                    "chunks": chunks,
                    "compressor": zarr.Blosc(),
                }
                kwargs["compute"] = False
                if overwrite:
                    kwargs["mode"] = "w"
        if kwargs["encoding"]:
            kwargs["encoding"]["time"] = {"dtype": "int32"}

    except KeyError:
        logging.error("Unable to encode chunks. Verify dataset.")
        raise

    return getattr(ds, f"to_{output_format}")(outfile, **kwargs)


def get_time_attrs(file_or_dataset: str | os.PathLike | xr.Dataset) -> (str, int):
    """Determine attributes related to time dimensions."""
    if isinstance(file_or_dataset, (str, Path)):
        ds = xr.open_dataset(Path(file_or_dataset).expanduser())
    else:
        ds = file_or_dataset

    calendar = ds.time.dt.calendar
    time = len(ds.time)

    return calendar, time


def get_global_attrs(
    file_or_dataset: str | os.PathLike | xr.Dataset,
) -> dict[str, str | int]:
    """Collect global attributes from NetCDF, Zarr, or Dataset object."""
    if isinstance(file_or_dataset, (str, Path)):
        file = Path(file_or_dataset).expanduser()
    elif isinstance(file_or_dataset, xr.Dataset):
        file = file_or_dataset
    else:
        raise NotImplementedError(f"Type: `{type(file_or_dataset)}`.")

    if isinstance(file, Path):
        if file.is_file() and file.suffix in [".nc", ".nc4"]:
            with nc.Dataset(file, mode="r") as ds:
                data = dict()
                for k in ds.ncattrs():
                    data[k] = getattr(ds, k)
        elif file.is_dir() and file.suffix == ".zarr":
            with zarr.open(file, mode="r") as ds:  # noqa
                data = ds.attrs.asdict()
    else:
        data = file.attrs

    return data


def sort_variables(
    files: list[Path], variables: Sequence[str]
) -> dict[str, list[Path]]:
    """Sort all variables within supplied files for treatment.

    Parameters
    ----------
    files: list of Path
    variables: sequence of str

    Returns
    -------
    dict[str, list[Path]]
    """
    variable_sorted = dict()
    if variables:
        logging.info("Sorting variables into groups. This could take some time.")
        for variable in variables:
            var_group = []
            for file in files:
                if file.name.startswith(variable):
                    var_group.append(file)
            if not var_group:
                logging.warning(f"No files found for {variable}. Continuing...")
                continue
            variable_sorted[variable] = sorted(var_group)
    else:
        variable_sorted["all_variables"] = files

    return variable_sorted


def get_chunks_on_disk(file: os.PathLike | str) -> dict:
    """Determine the chunks on disk for a given NetCDF or Zarr file.

    Parameters
    ----------
    file : str or os.PathLike
        File to be examined. Supports NetCDF and Zarr.

    Returns
    -------
    dict
    """
    chunks = dict()
    file = Path(file)

    if file.suffix.lower() in [".nc", ".nc4"]:
        with nc.Dataset(file) as ds:
            for v in ds.variables:
                chunks[v] = dict()
                for ii, dim in enumerate(ds[v].dimensions):
                    chunks[v][dim] = ds[v].chunking()[ii]
    elif file.suffix.lower() == "zarr" and file.is_dir():
        with zarr.open(file, "r") as ds:  # noqa
            for v in ds.arrays():
                # Check if variable is chunked
                if v[1]:
                    chunks[v[0]] = v[1]
    else:
        raise NotImplementedError(f"File type: {file.suffix}.")
    return chunks


def creation_date(path_to_file: str | os.PathLike) -> float | date:
    """Return the date that a file was created, falling back to when it was last modified if unable to determine.

    See https://stackoverflow.com/a/39501288/1709587 for explanation.

    Parameters
    ----------
    path_to_file : str or os.PathLike

    Returns
    -------
    float or date
    """
    if os.name == "nt":
        return Path(path_to_file).stat().st_ctime

    stat = Path(path_to_file).stat()
    try:
        return date.fromtimestamp(stat.st_ctime)
    except AttributeError:
        # We're probably on Linux. No easy way to get creation dates here,
        # so we'll settle for when its content was last modified.
        return date.fromtimestamp(stat.st_mtime)
