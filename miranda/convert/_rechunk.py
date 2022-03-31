import logging.config
import os
import shutil
import time
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import xarray as xr
import zarr

from miranda.scripting import LOGGING_CONFIG

from . import reanalysis_project_institutes
from ._data_definitions import era5_variables

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["rechunk_reanalysis"]


def _rechunk_configurator(project, time_step):
    # ~35 Mo chunks
    if project.lower() in ["era5-single-levels", "era5", "era5-land"]:
        if time_step == "1hr":
            target_chunks = {
                "time": 24 * 7,
                "latitude": 225,
                "longitude": 252,
            }
        elif time_step == "day":
            target_chunks = {"time": 365, "latitude": 125, "longitude": 125}
        else:
            raise NotImplementedError()
    else:
        raise NotImplementedError()
    return target_chunks


def rechunk_reanalysis(
    project: str,
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    time_step: Optional[str] = None,
    target_chunks: Optional[Dict[str, int]] = None,
    variables: Optional[Sequence[str]] = None,
    output_format: str = "zarr",
    overwrite: bool = False,
):
    """Rechunks ERA5 dataset for better loading/reading performance.

    Parameters
    ----------
    project : {"era5", "era5-land", "era5-single-levels"}
    input_folder : str or os.PathLike
    output_folder : str or os.PathLike
    time_step : {"1hr", "day"}, optional
      Time step of the input data. Parsed from filename if not set.
    target_chunks : dict
      Must include "time", optionally "latitude" and "longitude"
    variables : Sequence[str]
    output_format : {"netcdf", "zarr"}
      Default: "zarr".
    overwrite : bool

    Returns
    -------
    None
    """
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)

    if time_step is None:
        test_file = next(input_folder.glob("*"))
        file_parts = str(test_file.name).split("_")
        time_step = file_parts[1]
        if time_step not in ["1hr", "day"]:
            raise NotImplementedError()

    if project.startswith("era5") and variables is None:
        variables = era5_variables.copy()

    errored = list()
    start_all = time.perf_counter()
    for variable in variables:
        try:
            next(input_folder.glob(f"{variable}*"))
        except StopIteration:
            logging.warning("No files found for %s. Continuing..." % variable)
            continue

        # STEP 1 : Rewrite all years chunked on spatial dimensions, but not along the time dimension.
        start_var = time.perf_counter()

        for file in sorted(
            input_folder.glob(
                f"{variable}_{time_step}_{reanalysis_project_institutes[project]}_{project}_reanalysis_*.nc"
            )
        ):
            start = time.perf_counter()

            if output_format == "netcdf":
                output_folder.mkdir(exist_ok=True)
                out = output_folder / f"{file.stem}.nc"
            elif output_format == "zarr":
                output_path = output_folder / "temp"
                output_path.mkdir(exist_ok=True, parents=True)
                out = output_path / f"{file.stem}.zarr"
            else:
                raise NotImplementedError()

            if (out.is_dir() or out.is_file()) and not overwrite:
                logging.info(f"Already completed: {file.name}")
                continue
            if out.is_dir() and output_format == "zarr":
                logging.warning(f"Removing existing zarr files for {out.name}.")
                shutil.rmtree(out)

            ds = xr.open_dataset(
                file,
                chunks={"time": -1},
            )

            if target_chunks is None:
                target_chunks = _rechunk_configurator(project, time_step)

            # Set correct chunks in encoding options
            encoding = dict()
            try:
                for name, da in ds.data_vars.items():
                    chunks = list()
                    for dim in da.dims:
                        if dim in target_chunks.keys():
                            chunks.append(target_chunks[str(dim)])
                        else:
                            chunks.append(len(da[dim]))

                    if output_format == "netcdf":
                        encoding[name] = {
                            "chunksizes": chunks,
                            "zlib": True,
                        }
                    elif output_format == "zarr":
                        encoding[name] = {
                            "chunks": chunks,
                            "compressor": zarr.Blosc(),
                        }
            except KeyError:
                logging.warning(f"{file} has chunking errors. Verify data manually.")
                errored.append(file)
                continue

            # Write yearly file
            getattr(ds, f"to_{output_format}")(out, encoding=encoding)

            logging.info(f"Done for {file.stem} in {time.perf_counter() - start:.2f} s")

        if output_format == "netcdf":
            continue

        logging.info(
            f"First step done for {variable} in {(time.perf_counter() - start_all) / 3600:.2f} h"
        )

        # STEP 2 : Merge all years, chunking along the time dimension.
        start = time.perf_counter()

        files = sorted(
            (output_folder / "temp").glob(
                f"{variable}_{time_step}_{reanalysis_project_institutes[project]}_{project}_reanalysis_*.zarr"
            )
        )

        ds = xr.open_mfdataset(files, parallel=True, engine="zarr")

        if time_step == "1hr":
            # Four months of hours
            ds = ds.chunk(dict(time=2922))
        elif time_step == "day":
            # Five years of days
            ds = ds.chunk(dict(time=1825))
        else:
            raise NotImplementedError()

        for var in ds.data_vars.values():
            del var.encoding["chunks"]

        merged_zarr = Path(
            output_folder
            / f"{variable}_{time_step}_{reanalysis_project_institutes[project]}_{project}_reanalysis.zarr"
        )
        try:
            ds.to_zarr(merged_zarr, mode="w" if overwrite else "w-")
        except zarr.errors.ContainsGroupError:
            logging.error(
                'Files exist for variable %s. Consider using "overwrite=True"'
                % variable
            )
            raise

        logging.info(
            f"Second step done for {variable} in {time.perf_counter() - start:.2f} s"
        )
        logging.info(
            f"Both steps done for {variable} in {time.perf_counter() - start_var:.2f} s"
        )

    logging.info(
        f"All variables in {input_folder} done in {time.perf_counter() - start_all}"
    )
    if errored:
        errored_filenames = "\n".join(map(str, errored))
        logging.warning(f"Errored files are as follows: \n{errored_filenames}")
