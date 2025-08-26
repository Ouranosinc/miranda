"""IO Output Operations module."""

from __future__ import annotations
import logging
import os
import shutil
import time
import warnings
from collections.abc import Sequence
from pathlib import Path

import dask
import xarray as xr
from dask.distributed import Client

from miranda.convert.utils import date_parser

from ._input import discover_data
from ._rechunk import fetch_chunk_config, prepare_chunks_for_ds, translate_time_chunk
from .utils import delayed_write, get_global_attrs, name_output_file, sort_variables


logger = logging.getLogger("miranda.io.output")

__all__ = [
    "concat_rechunk_zarr",
    "merge_rechunk_zarrs",
    "write_dataset",
    "write_dataset_dict",
]


def write_dataset(
    ds: xr.DataArray | xr.Dataset,
    output_path: str | os.PathLike,
    output_format: str,
    output_name: str | None = None,
    chunks: dict | None = None,
    overwrite: bool = False,
    compute: bool = True,
) -> dict[str, Path]:
    """
    Write xarray object to NetCDf or Zarr with appropriate chunking regime.

    Parameters
    ----------
    ds : xr.DataArray or xr.Dataset
        Dataset or DatArray.
    output_path : str or os.PathLike
        Output folder path.
    output_format: {"netcdf", "zarr"}
        Output data container type.
    output_name: str, optional
        Output file name.
    chunks : dict, optional
        Chunking layout to be written to new files. If None, chunking will be left to the relevant backend engine.
    overwrite : bool
        Whether to remove existing files or fail if files already exist.
    compute : bool
        If True, files will be converted with each call to file conversion.
        If False, will return a dask.Delayed object that can be computed later.
        Default: True.

    Returns
    -------
    dict[str, Path]
    """
    if isinstance(output_path, str):
        output_path = Path(output_path)

    if not output_name:
        output_name = name_output_file(ds, output_format)
    else:
        output_name = str(output_name)

    outfile_path = output_path.joinpath(output_name)

    if overwrite and outfile_path.exists():
        msg = f"Removing existing {output_format} files for {outfile_path}."
        warnings.warn(msg, stacklevel=2)
        if outfile_path.is_dir():
            shutil.rmtree(outfile_path)
        if outfile_path.is_file():
            outfile_path.unlink()

    if chunks is None and "frequency" in ds.attrs:
        freq = ds.attrs.get("frequency")
        if not freq:
            raise ValueError("If 'chunks' are not provided, the 'frequency' attribute must be set.")
        if "lat" in ds.dims and "lon" in ds.dims:
            chunks = fetch_chunk_config(priority="time", freq=freq, dims=ds.dims)
        elif "lat" not in ds.dims and "lon" not in ds.dims:
            chunks = fetch_chunk_config(priority="stations", freq=freq, dims=ds.dims)

    msg = f"Writing {outfile_path}."
    logger.info(msg)
    write_object = delayed_write(
        ds,
        outfile_path,
        output_format,
        overwrite,
        target_chunks=prepare_chunks_for_ds(ds, chunks),
    )
    if compute:
        write_object.compute()
        return dict(path=outfile_path)
    return dict(path=outfile_path, object=write_object)


def write_dataset_dict(
    dataset_dict: dict[str, xr.Dataset | None],
    output_folder: str | os.PathLike,
    temp_folder: str | os.PathLike,
    *,
    output_format: str = "zarr",
    overwrite: bool = False,
    chunks: dict[str, int],
    **dask_kwargs,
):
    r"""
    Write dataset from Miranda-formatted dataset.

    Parameters
    ----------
    dataset_dict : dict[str, xr.Dataset or None]
    output_folder : str or os.PathLike
    temp_folder : str or os.PathLike
    output_format : {"netcdf", "zarr"}
    overwrite : bool
    chunks : dict[str, int]
    \*\*dask_kwargs

    Returns
    -------
    None
    """
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()
    if isinstance(temp_folder, str):
        temp_folder = Path(temp_folder).expanduser()

    if temp_folder:
        if temp_folder.exists():
            shutil.rmtree(temp_folder)
        temp_folder.mkdir(parents=True)

    for variable, ds in dataset_dict.items():
        if ds is None:
            continue

        outfile = name_output_file(ds, output_format)
        outpath = output_folder.joinpath(str(variable), outfile)
        if not outpath.exists() or overwrite:
            outpath.parent.mkdir(parents=True, exist_ok=True)
            if outpath.exists():
                shutil.rmtree(outpath)

            tmp_path = None
            if temp_folder:
                tmp_path = temp_folder.joinpath(outfile)
            job = delayed_write(
                ds,
                tmp_path if tmp_path else outpath,
                output_format=output_format,
                target_chunks=chunks,
                overwrite=True,
            )
            with Client(**dask_kwargs):
                dask.compute(job)
            if temp_folder:
                shutil.copytree(tmp_path, outpath, dirs_exist_ok=True)
                shutil.rmtree(tmp_path)

        elif outpath.exists() and not overwrite:
            existing_ds = xr.open_dataset(outpath, engine=output_format, decode_times=False)
            if "time" in ds.dims and "time" in existing_ds.dims:
                if ds.time.size > existing_ds.time.size:
                    msg = f"Dataset {variable} has more time points than existing file. Will overwrite {outpath.as_posix()}."
                    warnings.warn(msg, stacklevel=2)

                    tmp_path = None
                    if temp_folder:
                        tmp_path = temp_folder.joinpath(outfile)
                    job = delayed_write(
                        ds,
                        tmp_path if tmp_path else outpath,
                        output_format=output_format,
                        target_chunks=chunks,
                        overwrite=True,
                    )
                    with Client(**dask_kwargs):
                        dask.compute(job)
                    if temp_folder:
                        shutil.copytree(tmp_path, outpath, dirs_exist_ok=True)
                        shutil.rmtree(tmp_path)

        else:
            msg = f"Skipping {outpath.as_posix()} as overwrite is False and time dimension is sufficient."
            warnings.warn(msg, stacklevel=2)


# FIXME: concat_rechunk and merge_rechunk could be collapsed into each other
def concat_rechunk_zarr(
    freq: str,
    input_folder: str | os.PathLike,
    output_folder: str | os.PathLike,
    overwrite: bool = False,
    **dask_kwargs,
) -> None:
    r"""
    Concatenate and rechunk zarr files.

    Parameters
    ----------
    freq : str
    input_folder : str or os.PathLike
    output_folder : str or os.PathLike
    overwrite : bool
    \*\*dask_kwargs

    Returns
    -------
    None
    """
    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()

    list_zarr = sorted(list(input_folder.glob("*.zarr")))

    out_stem = "_".join(list_zarr[0].stem.split("_")[0:-1])
    start_year = date_parser(list_zarr[0].stem.split("_")[-1], output_type="datetime").year
    end_year = date_parser(list_zarr[-1].stem.split("_")[-1], output_type="datetime").year

    out_zarr = output_folder.joinpath(f"{out_stem}_{start_year}_{end_year}.zarr")

    if not out_zarr.exists() or overwrite:
        if out_zarr.exists():
            shutil.rmtree(out_zarr)
        # maketemp files 1 zarr per 4 years
        years = [y for y in range(int(start_year), int(end_year) + 1)]
        years = [years[x : x + 4] for x in range(0, len(years), 4)]
        tmp_folder = out_zarr.parent.joinpath("tmp")
        tmp_folder.mkdir(parents=True, exist_ok=True)
        chunks = dict()
        for year in years:
            list_zarr1 = sorted([zarr_file for zarr_file in list_zarr if int(zarr_file.stem.split("_")[-1].split("-")[0][0:4]) in year])
            # assert len(list_zarr1) / len(year) == 12
            ds = xr.open_mfdataset(list_zarr1, parallel=True, engine="zarr")
            if not chunks:
                chunk_config = fetch_chunk_config(priority="time", freq=freq, dims=ds.dims)
                chunks.update(
                    translate_time_chunk(
                        chunks=chunk_config,
                        calendar=ds.time.dt.calendar,
                        timesize=len(ds.time),
                    )
                )
            tmp_zarr = tmp_folder.joinpath(
                f"{out_zarr.stem.split(f'_{start_year}_')[0]}_{year[0]}-{year[-1]}.zarr",
            )
            msg = f"Writing year {year} to {tmp_zarr.as_posix()}."
            logger.info(msg)

            job = delayed_write(
                ds=ds,
                outfile=tmp_zarr,
                output_format="zarr",
                target_chunks=chunks,
                overwrite=True,
            )  # kwargs=zarr_kwargs)
            # FIXME: Client is only needed for computation. Should be elsewhere.
            with Client(**dask_kwargs):
                dask.compute(job)

        # write to final
        tmpzarrlist = sorted(list(tmp_zarr.parent.glob("*.zarr")))
        del ds
        ds = xr.open_mfdataset(tmpzarrlist, engine="zarr", parallel=True, combine="by_coords")
        zarr_kwargs = None
        job = delayed_write(
            ds=ds,
            outfile=out_zarr,
            output_format="zarr",
            target_chunks=chunks,
            overwrite=False,
            encode=False,
            kwargs=zarr_kwargs,
        )
        with Client(**dask_kwargs):
            dask.compute(job)
        shutil.rmtree(tmp_folder)


def merge_rechunk_zarrs(
    input_folder: str | os.PathLike,
    output_folder: str | os.PathLike,
    project: str | None = None,
    target_chunks: dict[str, int] | None = None,
    variables: Sequence[str] | None = None,
    freq: str | None = None,
    suffix: str = "zarr",
    overwrite: bool = False,
) -> None:
    """
    Merge and rechunk zarr files.

    Parameters
    ----------
    input_folder : str or os.PathLike
    output_folder : str or os.PathLike
    project : str, optional
    target_chunks : dict[str, int], optional
    variables : Sequence of str, optional
    freq : str, optional
    suffix : {"nc", "zarr"}
    overwrite : bool

    Returns
    -------
    None
    """
    chunk_defaults = {
        # Four months of hours
        "1hr": {"time": 2922},
        # Four years of days
        "day": {"time": 365 * 4 + 1},
    }

    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()
    output_folder.mkdir(exist_ok=True)

    files_found = discover_data(input_folder, suffix=suffix)
    variable_sorted = sort_variables(files_found, variables)

    if project is None and target_chunks is None:
        warnings.warn("`project` and `target_chunks` not set. Attempting to find `project` from attributes", stacklevel=2)
        project = get_global_attrs(files_found[0]).get("project")
        if not project:
            raise ValueError("`project` not found. Must pass either `project` or `target_chunks`.")

    if not freq:
        freq = get_global_attrs(files_found[0]).get("frequency")
        if not freq:
            raise ValueError("Frequency not found in file attributes.")
    if not target_chunks:
        try:
            target_chunks = chunk_defaults[freq]
        except KeyError as err:
            msg = f"Frequency not supported: `{freq}`."
            raise NotImplementedError(msg) from err

    start = time.perf_counter()

    for variable in variable_sorted.keys():
        start_var = time.perf_counter()

        ds = xr.open_mfdataset(files_found, parallel=True, engine="zarr")

        merged_zarr = name_output_file(ds, output_format="zarr")
        out_zarr = output_folder.joinpath(merged_zarr)

        if overwrite:
            if out_zarr.is_dir():
                msg = f"Removing existing zarr files for {out_zarr.name}."
                warnings.warn(msg, stacklevel=2)
                shutil.rmtree(out_zarr)
        else:
            msg = f"Files exist: {out_zarr.name}. Skipping..."
            logger.info(msg)
            continue

        ds = ds.chunk(target_chunks)
        for var in ds.data_vars.values():
            del var.encoding["chunks"]
        ds.to_zarr(out_zarr, mode="w")

        msg = f"{variable} rechunked in {(time.perf_counter() - start_var) / 3600:.2f} h."
        logger.info(msg)
    msg = f"All variables rechunked in {time.perf_counter() - start:.2f} s."
    logger.info(msg)
