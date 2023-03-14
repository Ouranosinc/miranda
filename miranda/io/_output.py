import logging.config
import os
import shutil
from pathlib import Path
from typing import Optional, Union

import xarray as xr
import zarr
from dask import compute, delayed
from dask.distributed import Client

from miranda.decode import date_parser
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["name_output_file", "delayed_write", "write_dataset", "concat_zarr"]


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
) -> delayed:
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


def concat_zarr(
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
        input_folder = Path(input_folder)
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

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
                compute(job)

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
            compute(job)

        shutil.rmtree(tmpzarr.parent)
