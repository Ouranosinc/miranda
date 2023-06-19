"""Specialized conversion tools for Environment and Climate Change Canada / Meteorological Service of Canada data."""

from __future__ import annotations

import json
import logging.config
import multiprocessing as mp
import os
import time
from functools import partial
from pathlib import Path

from miranda.eccc._raw import _convert_station_file
from miranda.eccc._utils import cf_station_metadata
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


_data_folder = Path(__file__).parent / "data"

eccc_observation_variables = dict()
eccc_observation_variables["flat"] = [
    v
    for v in json.load(open(_data_folder / "eccc_obs_flat_attrs.json"))[
        "variables"
    ].keys()
]
eccc_observation_variables["summary"] = [
    attrs["_cf_variable_name"]
    for attrs in json.load(open(_data_folder / "eccc_obs_summary_cf_attrs.json"))[
        "variables"
    ].values()
]
eccc_observation_variables["homogenized"] = [
    attrs["_cf_variable_name"]
    for attrs in json.load(open(_data_folder / "eccc_homogenized_cf_attrs.json"))[
        "variables"
    ].values()
]


def convert_flat_files(
    source_files: str | os.PathLike,
    output_folder: str | os.PathLike | list[str | int],
    variables: str | int | list[str | int],
    mode: str = "hourly",
    n_workers: int = 4,
) -> None:
    """

    Parameters
    ----------
    source_files: str or Path
    output_folder: str or Path
    variables: str or List[str]
    mode: {"hourly", "daily"}
    n_workers: int

    Returns
    -------
    None
    """
    func_time = time.time()

    if mode.lower() in ["h", "hour", "hourly"]:
        num_observations = 24
        column_names = ["code", "year", "month", "day", "code_var"]
        column_dtypes = [str, float, float, float, str]
    elif mode.lower() in ["d", "day", "daily"]:
        num_observations = 31
        column_names = ["code", "year", "month", "code_var"]
        column_dtypes = [str, float, float, str]
    else:
        raise NotImplementedError("`mode` must be 'h'/'hourly or 'd'/'daily'.")

    # Preparing the data column headers
    for i in range(1, num_observations + 1):
        data_entry, flag_entry = f"D{i:0n}", f"F{i:0n}"
        column_names.append(data_entry)
        column_names.append(flag_entry)
        column_dtypes.extend([str, str])

    if isinstance(variables, (str, int)):
        variables = [variables]

    for variable_code in variables:
        variable_code = str(variable_code).zfill(3)
        metadata = cf_station_metadata(variable_code)
        nc_name = metadata["nc_name"]

        rep_nc = Path(output_folder).joinpath(nc_name)
        rep_nc.mkdir(parents=True, exist_ok=True)

        # Loop on the files
        logging.info(
            f"Collecting files for variable '{metadata['standard_name']}' "
            f"(filenames containing '{metadata['_table_name']}')."
        )
        list_files = list()
        if isinstance(source_files, list) or Path(source_files).is_file():
            list_files.append(source_files)
        else:
            glob_patterns = [g for g in metadata["_table_name"]]
            for pattern in glob_patterns:
                list_files.extend(
                    [f for f in Path(source_files).rglob(f"{pattern}*") if f.is_file()]
                )
        manager = mp.Manager()
        errored_files = manager.list()
        converter_func = partial(
            _convert_station_file,
            output_path=rep_nc,
            errored_files=errored_files,
            mode=mode,
            variable_code=variable_code,
            column_names=column_names,
            column_dtypes=column_dtypes,
            **metadata,
        )
        with mp.Pool(processes=n_workers) as pool:
            pool.map(converter_func, list_files)
            pool.close()
            pool.join()

        if errored_files:
            logging.warning(
                "Some files failed to be properly parsed:\n", ", ".join(errored_files)
            )

    logging.warning(f"Process completed in {time.time() - func_time:.2f} seconds")
