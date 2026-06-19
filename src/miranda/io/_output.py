"""IO Output Operations module."""

from __future__ import annotations
import logging
import os
import shutil
import warnings
from pathlib import Path
from typing import Literal

import dask
import xarray as xr
from dask.diagnostics import ProgressBar
from dask.distributed import Client

from ._rechunk import fetch_chunk_config, prepare_chunks_for_ds
from .utils import delayed_write, name_output_file


logger = logging.getLogger("miranda.io.output")

__all__ = ["write_dataset", "write_dataset_dict", "write_zarr"]


def write_dataset(
    ds: xr.DataArray | xr.Dataset,
    output_path: str | os.PathLike,
    output_format: Literal["netcdf", "zarr"],
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
    dict[str, Path | dask.delayed.Delayed]
        Dictionary with either {"path": outfile_path} or {"path": path, "object": write_object} for job writing.
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
    """
    Write dataset from Miranda-formatted dataset.

    Parameters
    ----------
    dataset_dict : dict[str, xr.Dataset or None]
    output_folder : str or os.PathLike
    temp_folder : str or os.PathLike
    output_format : {"netcdf", "zarr"}
    overwrite : bool
    chunks : dict[str, int]
    **dask_kwargs

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


def write_zarr(
    ds: xr.Dataset,
    out_zarr: Path,
    chunks: dict,
    zarr_format: int = 2,
    overwrite: bool = False,
) -> None:
    """
    Write Zarr file.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset.
    out_zarr : Path
        Output Zarr file.
    chunks : dict
        Chunk sizes.
        If an empty dictionary is passed, no chunking will be performed.
    zarr_format : int
        Zarr format version (2 or 3). Default is 2.
    overwrite : bool
        Whether to overwrite. Default is False.
    """
    if not out_zarr.exists() or overwrite:
        with ProgressBar():
            for vv in ds.data_vars:
                if ds[vv].dtype == object:
                    ds[vv] = ds[vv].astype(str)
            if len(chunks):
                ds.chunk(chunks).to_zarr(out_zarr.with_suffix(".tmp.zarr"), mode="w", zarr_format=zarr_format)
            else:
                ds.to_zarr(out_zarr.with_suffix(".tmp.zarr"), mode="w", zarr_format=zarr_format)
        shutil.move(out_zarr.with_suffix(".tmp.zarr"), out_zarr)
