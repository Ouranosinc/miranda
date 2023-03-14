import logging.config
import os
import shutil
import time
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import xarray as xr
import zarr

from miranda.convert import project_institutes
from miranda.io import delayed_write
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["config_chunks", "rechunk_files"]

_chunk_presets = dict()
_chunk_presets[
    "era5-single-levels",
    "era5-single-levels-preliminary-back-extension",
    "era5-land",
    "era5-pressure-levels",
    "era5-pressure-levels-preliminary-back-extension",
] = {
    "1hr": {
        "time": 24 * 7,
        "latitude": 225,
        "longitude": 252,
    },
}


def config_chunks(project, time_step: str, levels: int = None):
    # ~35 Mo chunks
    if project.lower() in [
        "era5-single-levels",
        "era5-single-levels-preliminary-back-extension",
        "era5-land",
    ] or (project.lower() in [] and levels == 1):
        if time_step == "1hr":
            # Chunks for monthly files, optimized for Zarr
            target_chunks = {
                "time": 24 * 7,
                "latitude": 225,
                "longitude": 252,
            }
        elif time_step == "day":
            # Chunks for annual files, optimized for Zarr
            target_chunks = {"time": 365, "latitude": 125, "longitude": 125}
        else:
            raise NotImplementedError()
    else:
        raise NotImplementedError()
    return target_chunks


def rechunk_files(
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    project: str,
    time_step: Optional[str] = None,
    target_chunks: Optional[Dict[str, int]] = None,
    variables: Optional[Sequence[str]] = None,
    output_format: str = "zarr",
    overwrite: bool = False,
) -> None:
    """Rechunks ERA5 dataset for better loading/reading performance.

    Warnings
    --------
    Globbing assumes that target datasets to be rechunked have been saved in NetCDF format.
    File naming requires the following order of facets: `{variable}_{time_step}_{institute}_{project}_reanalysis_*.nc`.

    Parameters
    ----------
    input_folder : str or os.PathLike
      Folder to be examined. Performs globbing.
    output_folder : str or os.PathLike
      Target folder.
    project : {"era5-land", "era5-single-levels", "era5-single-levels-preliminary-back-extension", "era5-single-levels", "era5-pressure-levels-preliminary-back-extension"}
      Supported projects.
    time_step : {"1hr", "day"}, optional
      Time step of the input data. Parsed from filename if not set.
    target_chunks : dict
      Must include "time", optionally "latitude" and "longitude"
    variables : Sequence[str], optional
      If no variables set, will attempt to process all variables supported based on project name.
    output_format : {"netcdf", "zarr"}
      Default: "zarr".
    overwrite : bool
      Will overwrite files. For zarr, existing folders will be removed before writing.

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
            raise NotImplementedError(f"Time step: `{time_step}`.")

    errored = list()
    start_all = time.perf_counter()

    glob_patterns = []
    if variables:
        for variable in variables:
            try:
                next(input_folder.glob(f"{variable}*"))
            except StopIteration:
                logging.warning(f"No files found for {variable}. Continuing...")
                continue
            glob_patterns.append(
                f"{variable}_{time_step}_{project_institutes[project]}_{project}_*.nc"
            )
    else:
        glob_patterns.append(
            f"*_{time_step}_{project_institutes[project]}_{project}_*.nc"
        )

    for pattern in glob_patterns:
        start_var = time.perf_counter()
        variable = pattern.split("_")[0]

        # STEP 1 : Rewrite all years chunked on spatial dimensions, but not along the time dimension.
        for file in sorted(input_folder.glob(pattern)):
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

            ds = xr.open_dataset(file, chunks={"time": -1})

            try:
                levels = len(ds.level)
            except AttributeError:
                levels = None

            if target_chunks is None:
                target_chunks = config_chunks(project, time_step, levels)

            try:
                job = delayed_write(
                    ds,
                    out,
                    output_format=output_format,
                    overwrite=overwrite,
                    target_chunks=target_chunks,
                )

            except KeyError:
                logging.warning(f"{file} has chunking errors. Verify data manually.")
                errored.append(file)
                continue

            # Write out rechunked files at the same output frequencies as they were read.
            job.compute()
            logging.info(f"Done for {file.stem} in {time.perf_counter() - start:.2f} s")

        if output_format == "netcdf":
            logging.info(
                f"{variable} rechunked in {(time.perf_counter() - start_all) / 3600:.2f} h"
            )
            continue

        logging.info(
            f"First step done for {variable} in {(time.perf_counter() - start_all) / 3600:.2f} h"
        )

        # STEP 2 : Merge all years, chunking along the time dimension.
        start = time.perf_counter()

        files = sorted(
            (output_folder / "temp").glob(
                f"{variable}_{time_step}_{project_institutes[project]}_{project}_*.zarr"
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
            / f"{variable}_{time_step}_{project_institutes[project]}_{project}.zarr"
        )
        try:
            ds.to_zarr(merged_zarr, mode="w" if overwrite else "w-")
        except zarr.errors.ContainsGroupError:
            logging.error(
                f"Files exist for pattern {pattern.split('_')[0]}. Consider using `overwrite=True`"
            )
            raise

        logging.info(
            f"Second step done for {pattern.split('_')[0]} in {time.perf_counter() - start:.2f} s"
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
