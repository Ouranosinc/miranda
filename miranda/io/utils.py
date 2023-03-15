import logging.config
import os
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import dask
import netCDF4 as nc  # noqa
import xarray as xr
import zarr

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = [
    "creation_date",
    "delayed_write",
    "get_attributes",
    "get_chunks_on_disk",
    "name_output_file",
    "sort_variables",
]


def name_output_file(ds: xr.Dataset, project: str, output_format: str) -> str:
    """

    Parameters
    ----------
    ds : xr.Dataset
    project: str
    output_format : {"netcdf", "zarr"}
        Suffix to be used for filename

    Returns
    -------
    str
    """
    if output_format.lower() not in {"netcdf", "zarr"}:
        raise NotImplementedError(f"Format: {output_format}.")
    else:
        suffix = dict(netcdf="nc", zarr="zarr")[output_format]

    var_name = list(ds.data_vars.keys())[0]
    time_freq = ds.attrs.get("frequency")
    institution = ds.attrs.get("institution")
    time_start, time_end = ds.time.isel(time=[0, -1]).dt.strftime("%Y%m%d").values

    return f"{var_name}_{time_freq}_{institution}_{project}_{time_start}-{time_end}.{suffix}"


def delayed_write(
    ds: xr.Dataset,
    outfile: Union[str, os.PathLike],
    output_format: str,
    overwrite: bool,
    target_chunks: Optional[dict] = None,
) -> dask.delayed:
    """

    Parameters
    ----------
    ds : Union[xr.Dataset, str, os.PathLike]
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


def get_attributes(file: Union[str, os.PathLike]):
    if isinstance(file, str):
        file = Path(file).expanduser()

    if file.is_file() and file.suffix in [".nc", ".nc4"]:
        with nc.Dataset(file, mode="r") as ds:
            data = dict()
            for k in ds.ncattrs():
                data[k] = getattr(ds, k)
    elif file.is_dir() and file.suffix == ".zarr":
        with zarr.open(file, mode="r") as ds:  # noqa
            data = ds.attrs.asdict()

    return data


def sort_variables(
    files: List[Path], variables: Sequence[str]
) -> Dict[str, List[Path]]:
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


def get_chunks_on_disk(file: Union[os.PathLike, str]) -> dict:
    """

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


def creation_date(path_to_file: Union[Path, str]) -> Union[float, date]:
    """Try to get the date that a file was created, falling back to when it was last modified if that isn't possible.

    See https://stackoverflow.com/a/39501288/1709587 for explanation.

    Parameters
    ----------
    path_to_file: Union[Path, str]

    Returns
    -------
    Union[float, date]
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
