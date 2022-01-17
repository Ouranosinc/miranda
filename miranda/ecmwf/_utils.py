import logging.config
import os
import time
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import xarray as xr
import zarr

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["rechunk_ecmwf"]


def rechunk_ecmwf(
    project: str,
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    time_step: str,
    target_chunks: Optional[Dict[str, int]] = None,
    variables: Optional[Sequence[str]] = None,
    output_format: str = ".nc",
):
    """Rechunks ERA5 dataset for better loading/reading performance.

    Parameters
    ----------
    project: {"era5", "era5-land", "era5-single-levels"}
    input_folder: str or os.PathLike
    output_folder: str or os.PathLike
    time_step: {"hourly", "daily"}
    target_chunks: dict
      Must include "time", optionally "latitude" and "longitude"
    variables: Sequence[str]
    output_format: {".nc", ".zarr"}

    Returns
    -------
    None
    """
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)

    if variables is None:
        variables = ["potevap", "pr", "prsn", "snd", "swe", "tas", "td", "uas", "vas"]

    if target_chunks is None:
        # ~35 Mo chunks
        if project.lower() == ["era5-single-levels", "era5", "era5-land"]:
            if time_step == "hourly":
                target_chunks = {"time": 24 * 7, "latitude": 225, "longitude": 252}
            elif time_step == "daily":
                target_chunks = {"time": 365, "latitude": 125, "longitude": 125}
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()

    errored = list()
    for folder, variables in (input_folder, variables):
        start_all = time.perf_counter()
        for variable in variables:

            # STEP 1 : Rewrite all years chunked on spatial dimensions, but not along time.
            start_var = time.perf_counter()
            for file in sorted(
                folder.glob(f"{variable}_{project}_reanalysis_{time_step}_*.nc")
            ):
                start = time.perf_counter()

                if output_format == "nc":
                    output_folder.mkdir(exist_ok=True)
                    out = output_folder / f"{file.stem}.nc"
                elif output_format == "zarr":
                    outpath = output_folder / "temp"
                    outpath.mkdir(exist_ok=True, parents=True)
                    out = outpath / f"{file.stem}.zarr"
                else:
                    raise NotImplementedError()

                if out.is_dir() or out.is_file():
                    logging.info(f"Already completed: {file}")
                    continue

                ds = xr.open_dataset(
                    file,
                    chunks={"time": -1},
                )

                # Set correct chunks in encoding options
                encoding = dict()
                try:
                    for name, da in ds.data_vars.items():
                        if output_format == "nc":
                            encoding[name] = {
                                "chunksizes": [
                                    target_chunks[str(dim)] for dim in da.dims
                                ],
                                "zlib": True,
                            }
                        elif output_format == "zarr":
                            encoding[name] = {
                                "chunks": [target_chunks[str(dim)] for dim in da.dims],
                                "compressor": zarr.Blosc(),
                            }
                except KeyError:
                    logging.warning(
                        f"{file} has chunking errors. Verify data manually."
                    )
                    errored.append(file)
                    continue

                # Write yearly file
                if output_format == "nc":
                    ds.to_netcdf(out, encoding=encoding)
                elif output_format == "zarr":
                    ds.to_zarr(out, encoding=encoding)

                logging.info(
                    f"Done for {file.stem} in {time.perf_counter() - start:.2f} s"
                )
            logging.info(
                f"First step done for {variable} in {(time.perf_counter() - start_all) / 3600:.2f} h"
            )

            if output_format == "nc":
                continue

            # STEP 2 : Merge all years, chunking along time.
            start = time.perf_counter()

            files = sorted(
                (output_folder / "temp").glob(
                    f"{variable}_{project}_reanalysis_{time_step}_*.zarr"
                )
            )

            ds = xr.open_mfdataset(files, parallel=True, engine="zarr")

            ds = ds.chunk({"time": 2922})
            for var in ds.data_vars.values():
                del var.encoding["chunks"]

            ds.to_zarr(
                output_folder / f"{variable}_{project}_reanalysis_{time_step}.zarr"
            )
            logging.info(
                f"Second step done for {variable} in {time.perf_counter() - start:.2f} s"
            )
            logging.info(
                f"Both steps done for {variable} in {time.perf_counter() - start_var:.2f} s"
            )

        print(f"All variables in {folder} done in {time.perf_counter() - start_all}")
        if errored:
            errored_filenames = "\n".join(map(str, errored))
            logging.warning(f"Errored files are as follows: \n{errored_filenames}")


def daily_metadata():
    pass


def hourly_metadata():
    pass


def correct_cf_metadata():
    pass


def increment_hourly_precipitation(
    files: Sequence[Union[os.PathLike, str]], output_folder
):
    pass
