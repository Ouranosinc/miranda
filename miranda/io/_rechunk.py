import json
import logging.config
import os
import shutil
import time
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import xarray as xr

from miranda.convert import project_institutes
from miranda.scripting import LOGGING_CONFIG

from ._input import discover_data
from .utils import delayed_write, get_global_attrs, get_time_attrs, sort_variables

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["fetch_chunk_config", "rechunk_files", "translate_time_chunk"]

_data_folder = Path(__file__).parent / "data"
chunks_configurations = json.load(open(_data_folder / "chunks.json"))


# Shamelessly stolen and modified from xscen
def translate_time_chunk(chunks: dict, calendar: str, timesize: int) -> dict:
    """Translate chunk specification for time into a number.

    Notes
    -----
    -1 translates to `timesize`
    'Nyear' translates to N times the number of days in a year of calendar `calendar`.
    """
    for k, v in chunks.items():
        if isinstance(v, dict):
            chunks[k] = translate_time_chunk(v.copy(), calendar, timesize)
        elif k == "time" and v is not None:
            if isinstance(v, str) and v.endswith("years"):
                n = int(str(chunks["time"]).strip("years").strip())
                Nt = n * {"noleap": 365, "360_day": 360, "all_leap": 366}.get(
                    calendar, 365.25
                )
                chunks[k] = int(Nt)
            elif v == -1:
                chunks[k] = timesize
    return chunks


def fetch_chunk_config(project: str, freq: str) -> Dict[str, int]:
    institute = project_institutes[project]
    entry = chunks_configurations[institute.upper()]

    # TODO: Currently no explicit handling for multi-level data
    if project in entry["projects"]:
        if freq in entry["chunks"]:
            return entry["chunks"][freq]
        raise KeyError(
            f"Chunks at frequency `{freq}` not found for project `{project}`."
        )
    raise KeyError(f"Project `{project}` not found.")


def rechunk_files(
    input_folder: Union[str, os.PathLike],
    output_folder: Union[str, os.PathLike],
    project: Optional[str] = None,
    time_step: Optional[str] = None,
    target_chunks: Optional[Dict[str, int]] = None,
    variables: Optional[Sequence[str]] = None,
    suffix: str = "nc",
    output_format: str = "zarr",
    overwrite: bool = False,
) -> None:
    """Rechunks dataset for better loading/reading performance.

    Warnings
    --------
    Globbing assumes that target datasets to be rechunked have been saved in NetCDF format.
    File naming requires the following order of facets: `{variable}_{time_step}_{institute}_{project}_reanalysis_*.nc`.
    Chunking dimensions are assumed to be CF-Compliant (`lat`, `lon`, `rlat`, `rlon`, `time`).

    Parameters
    ----------
    input_folder : str or os.PathLike
        Folder to be examined. Performs globbing.
    output_folder : str or os.PathLike
        Target folder.
    project : str, optional
        Supported projects. Used for determining chunk dictionary. Superseded if `target_chunks` is set.
    time_step : {"1hr", "day"}, optional
        Time step of the input data. Parsed from dataset attrs if not set. Superseded if `target_chunks` is set.
    target_chunks : dict, optional
        Must include "time", optionally "lat" and "lon", depending on dataset structure.
    variables : Sequence[str], optional
        If no variables set, will attempt to process all variables supported based on project name.
    suffix : {"nc", "zarr"}
        Suffix used to identify data files. Default: "nc".
    output_format : {"netcdf", "zarr"}
        Default: "zarr".
    overwrite : bool
        Will overwrite files. For zarr, existing folders will be removed before writing.

    Returns
    -------
    None
    """
    if isinstance(input_folder, str):
        input_folder = Path(input_folder).expanduser()
    if isinstance(output_folder, str):
        output_folder = Path(output_folder).expanduser()
    output_folder.mkdir(exist_ok=True)

    if project is None and target_chunks is None:
        logging.warning(
            "`project` and `target_chunks` not set. Attempting to find `project` from attributes"
        )
        project = get_global_attrs(next(input_folder.glob(f"*.{suffix}"))).get(
            "project"
        )
        if not project:
            raise ValueError(
                "`project` not found. Must pass either `project` or `target_chunks`."
            )

    if target_chunks is None:
        test_file = next(input_folder.glob(f"*.{suffix}"))
        if time_step is None:
            time_step = get_global_attrs(test_file).get("frequency")
            if not time_step:
                raise ValueError("Frequency not found in file attributes.")

    if output_format == "netcdf":
        output_suffix = "nc"
    elif output_format == "zarr":
        output_suffix = "zarr"
    else:
        raise NotImplementedError(f"Output format: `{output_format}`.")

    files_found = discover_data(input_folder, suffix=suffix)
    variable_sorted = sort_variables(files_found, variables)

    errored = list()
    start_all = time.perf_counter()
    for variable, files in variable_sorted.items():
        start_var = time.perf_counter()

        for file in files:
            start = time.perf_counter()

            output_path = output_folder / "temp"
            output_path.mkdir(exist_ok=True)
            out = output_path / f"{file.stem}.{output_suffix}"

            if out.is_dir() or out.is_file():
                if overwrite:
                    if out.is_dir() and output_format == "zarr":
                        logging.warning(f"Removing existing zarr files for {out.name}.")
                        shutil.rmtree(out)
                else:
                    logging.info(f"Files exist: {file.name}")
                    continue

            ds = xr.open_dataset(file, chunks={"time": -1})

            if target_chunks is None:
                chunk_config = fetch_chunk_config(project, time_step)
                calendar, size = get_time_attrs(ds)
                target_chunks = translate_time_chunk(chunk_config, calendar, size)

            try:
                delayed_write(
                    ds,
                    out,
                    output_format=output_format,
                    overwrite=overwrite,
                    target_chunks=target_chunks,
                ).compute()

            except KeyError:
                logging.warning(f"{file} has chunking errors. Verify data manually.")
                errored.append(file)
                continue

            logging.info(f"Done for {file.stem} in {time.perf_counter() - start:.2f} s")

            logging.info(f"Moving {file.name} to {output_folder}.")
            shutil.move(out, output_folder)

        logging.info(
            f"{variable} rechunked in {(time.perf_counter() - start_var) / 3600:.2f} h"
        )
        continue

    logging.info(
        f"All variables in {input_folder} rechunked in {time.perf_counter() - start_all}"
    )
    if errored:
        errored_filenames = "\n".join(map(str, errored))
        logging.warning(f"Errored files are as follows: \n{errored_filenames}")
