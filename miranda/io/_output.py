import logging.config
import os
import shutil
import time
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import dask
import xarray as xr
from dask.distributed import Client

from miranda.convert.utils import date_parser
from miranda.scripting import LOGGING_CONFIG

from ._input import discover_data
from .utils import delayed_write, get_attributes, name_output_file, sort_variables

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = [
    "write_dataset",
    "concat_rechunk_zarr",
    "merge_rechunk_zarrs",
]


def write_dataset(
    ds: Union[xr.DataArray, xr.Dataset],
    project: str,
    output_path: Union[str, os.PathLike],
    output_format: str,
    chunks: Optional[dict] = None,
    overwrite: bool = False,
    compute: bool = True,
):
    """

    Parameters
    ----------
    ds : xr.DataArray or xr.Dataset

    project : {"cordex", "cmip5", "cmip6", "ets-grnch", "isimip-ft", "pcic-candcs-u6", "converted"}
        Project name for decoding/handling purposes.
    output_path : str or os.PathLike
        Output folder path.
    output_format: {"netcdf", "zarr"}
        Output data container type.
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

    """
    if isinstance(output_path, str):
        output_path = Path(output_path)

    outfile = name_output_file(ds, project, output_format)
    outfile_path = output_path.joinpath(outfile)

    if overwrite and outfile_path.exists():
        logging.warning(f"Removing existing {output_format} files for {outfile}.")
        if outfile_path.is_dir():
            shutil.rmtree(outfile_path)
        if outfile_path.is_file():
            outfile_path.unlink()

    logging.info(f"Writing {outfile}.")
    write_object = delayed_write(
        ds,
        outfile_path,
        output_format,
        overwrite,
        target_chunks=chunks,
    )
    if compute:
        write_object.compute()
        return dict(path=outfile_path)
    return dict(path=outfile_path, object=write_object)


def write_dataset_dict(
    dataset_dict: Dict[str, xr.Dataset],
    output_folder: Union[str, os.PathLike],
    temp_folder: Union[str, os.PathLike],
    output_format: str = "zarr",
    overwrite: bool = False,
    chunks: Dict[str, int] = None,
    **dask_kwargs,
):
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

        outfile = name_output_file(ds, ds.attrs["project"], output_format)
        outpath = output_folder.joinpath(variable, outfile)
        if not outpath.exists() or overwrite:
            outpath.parent.mkdir(parents=True, exist_ok=True)
            if outpath.exists():
                shutil.rmtree(outfile)
            if temp_folder:
                tmp_path = temp_folder.joinpath(outfile)
                job = delayed_write(
                    ds,
                    tmp_path,
                    output_format=output_format,
                    target_chunks=chunks,
                    overwrite=True,
                )
                with Client(**dask_kwargs):
                    dask.compute(job)
                shutil.copytree(tmp_path, outpath, dirs_exist_ok=True)
                shutil.rmtree(tmp_path)
            else:
                job = delayed_write(
                    ds,
                    outpath,
                    output_format=output_format,
                    target_chunks=chunks,
                    overwrite=True,
                )
                with Client(**dask_kwargs):
                    dask.compute(job)
        else:
            logging.warning(
                f"{outpath.as_posix()} exists and overwrite is False. Continuing..."
            )


# FIXME: concat_rechunk and merge_rechunk could be collapsed into each other
def concat_rechunk_zarr(
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    overwrite: bool = False,
    **dask_kwargs,
):
    """

    Parameters
    ----------
    input_folder
    output_folder
    overwrite
    dask_kwargs

    Returns
    -------

    """
    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()

    list_zarr = sorted(list(input_folder.glob("*.zarr")))

    out_stem = "_".join(list_zarr[0].stem.split("_")[0:-1])
    start_year = date_parser(
        list_zarr[0].stem.split("_")[-1], output_type="datetime"
    ).year
    end_year = date_parser(
        list_zarr[-1].stem.split("_")[-1], output_type="datetime"
    ).year

    outzarr = f"{out_stem}_{start_year}_{end_year}.zarr"
    outzarr = output_folder.joinpath(outzarr)

    if not outzarr.exists() or overwrite:
        if "day" in input_folder.as_posix():
            chunks = dict(time=(365 * 4) + 1, rlon=50, rlat=50)
        else:
            chunks = dict(time=(24 * 30 * 2), rlon=50, rlat=50)

        # maketemp files 1 zarr per 4 years
        years = [y for y in range(int(start_year), int(end_year) + 1)]
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

            # FIXME: Client is only needed for computation. Should be elsewhere.
            with Client(**dask_kwargs):
                # if outzarr.exists():
                #     zarr_kwargs = {"append_dim": "time", "consolidated": True}
                # else:
                #     zarr_kwargs = {"consolidated": True}
                tmpzarr = outzarr.parent.joinpath(
                    "tmp",
                    f"{outzarr.stem.split(f'_{start_year}_')[0]}_{year[0]}-{year[-1]}.zarr",
                )
                tmpzarr.parent.mkdir(exist_ok=True, parents=True)
                logging.info(f"Writing year {year} to {tmpzarr.as_posix()}.")

                job = delayed_write(
                    ds=ds,
                    outfile=tmpzarr,
                    output_format="zarr",
                    target_chunks=chunks,
                    overwrite=overwrite,
                )  # kwargs=zarr_kwargs)
                dask.compute(job)

        # get tmp zarrs
        list_zarr = sorted(list(tmpzarr.parent.glob("*zarr")))
        ds = xr.open_mfdataset(list_zarr, engine="zarr")
        # FIXME: Client is only needed for computation. Should be elsewhere.
        with Client(**dask_kwargs):
            job = delayed_write(
                ds=ds,
                outfile=outzarr,
                output_format="zarr",
                target_chunks=chunks,
                overwrite=overwrite,
            )  # kwargs=zarr_kwargs)
            dask.compute(job)

        shutil.rmtree(tmpzarr.parent)


def merge_rechunk_zarrs(
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    project: Optional[str] = None,
    target_chunks: Optional[Dict[str, int]] = None,
    variables: Optional[Sequence[str]] = None,
    freq: Optional[str] = None,
    suffix: str = "zarr",
    overwrite: bool = False,
):
    """

    Parameters
    ----------
    input_folder
    output_folder
    project
    target_chunks
    variables
    freq
    suffix
    overwrite

    Returns
    -------

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
        logging.warning(
            "`project` and `target_chunks` not set. Attempting to find `project` from attributes"
        )
        project = get_attributes(next(input_folder.glob(f"*.{suffix}"))).get("project")
        if not project:
            raise ValueError(
                "`project` not found. Must pass either `project` or `target_chunks`."
            )

    if not freq:
        freq = get_attributes(files_found[0]).get("frequency")
        if not freq:
            raise ValueError("Frequency not found in file attributes.")
    if not target_chunks:
        try:
            target_chunks = chunk_defaults[freq]
        except KeyError:
            raise NotImplementedError(f"Frequency not supported: `{freq}`.")

    start = time.perf_counter()

    for variable, files in variable_sorted.items():
        start_var = time.perf_counter()

        ds = xr.open_mfdataset(files_found, parallel=True, engine="zarr")

        merged_zarr = name_output_file(ds, project, output_format="zarr")
        out_zarr = output_folder.joinpath(merged_zarr)

        if overwrite:
            if out_zarr.is_dir():
                logging.warning(f"Removing existing zarr files for {out_zarr.name}.")
                shutil.rmtree(out_zarr)
        else:
            logging.info(f"Files exist: {out_zarr.name}")
            continue

        ds = ds.chunk(target_chunks)
        for var in ds.data_vars.values():
            del var.encoding["chunks"]
        ds.to_zarr(out_zarr, mode="w")

        logging.info(
            f"{variable} rechunked in {(time.perf_counter() - start_var) / 3600:.2f} h"
        )
    logging.info(f"All variables rechunked in {time.perf_counter() - start:.2f} s")
