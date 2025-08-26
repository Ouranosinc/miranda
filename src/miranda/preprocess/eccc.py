"""Specialized conversion tools for Environment and Climate Change Canada / Meteorological Service of Canada data."""

from __future__ import annotations
import contextlib
import logging
import tempfile
from collections.abc import Callable
from pathlib import Path

from dask.diagnostics import ProgressBar

from miranda.storage import file_size, report_file_size
from miranda.utils import generic_extract_archive


_data_folder = Path(__file__).parent / "configs"


def _run_func_on_archive_with_optional_dask(
    file: Path,
    function: Callable,
    errored_files: list[Path],
    **dask_kwargs,
) -> None:
    r"""
    Run a function on a file archive, extracting it if necessary.

    Parameters
    ----------
    file : Path
        File archive to process.
    function : Callable
        Function to run on the file.
    errored_files : list[Path]
        List of files that errored during processing.
    \*\*dask_kwargs : Any
        Keyword arguments to pass to dask.distributed.Client.

    Notes
    -----
    If the file is larger than 1 GiB or dask_kwargs are passed, dask.dataframes will be used.
    Partial function requires the function to accept the following parameters:
    - file: Path
    - using_dask: bool
    - client: dask.distributed.Client
    """
    with tempfile.TemporaryDirectory() as temp_folder:
        if file.suffix in [".gz", ".tar", ".zip", ".7z"]:
            data_files = generic_extract_archive(file, output_dir=temp_folder)
        else:
            data_files = [file]
        msg = f"Processing file: {file}."
        logging.info(msg)

        # 1 GiB
        size_limit = 2**30

        for data in data_files:
            size = file_size(data)
            if size > size_limit or dask_kwargs:
                if dask_kwargs:
                    logging.info("`dask_kwargs` provided - Using dask.dataframes.")
                elif size > size_limit:
                    msg = f"File exceeds {report_file_size(size_limit)} - Using dask.dataframes."
                    logging.info(msg)
                client = ProgressBar
                using_dask = True
            else:
                msg = f"File below {report_file_size(size_limit)} - Using pandas.dataframes."
                logging.info(msg)
                client = contextlib.nullcontext
                using_dask = False

            with client(**dask_kwargs) as c:
                try:
                    function(data, using_dask=using_dask, client=c)
                except FileNotFoundError:
                    errored_files.append(data)

        if Path(temp_folder).iterdir():
            for temporary_file in Path(temp_folder).glob("*"):
                if temporary_file in data_files:
                    temporary_file.unlink()


# def convert_flat_files(
#     source_files: str | os.PathLike,
#     output_folder: str | os.PathLike | list[str | int],
#     variables: str | int | list[str | int],
#     project: str = "eccc-obs",
#     mode: str = "hourly",
#     **dask_kwargs,
# ) -> None:
#     """
#
#     Parameters
#     ----------
#     source_files: str or Path
#     output_folder: str or Path
#     variables: str or List[str]
#     project: {"eccc-obs", "eccc-obs-summary", "eccc-homogenized"}
#     mode: {"hourly", "daily"}
#
#     Returns
#     -------
#     None
#     """
#
#     if isinstance(variables, (str, int)):
#         variables = [variables]
#
#     for variable_code in variables:
#         variable_code = str(variable_code).zfill(3)
#         metadata = load_json_data_mappings("eccc-obs").get(variable_code)
#
#
#
#         # Loop on the files
#         logging.info(
#             f"Collecting files for variable '{metadata['standard_name']}' "
#             f"(filenames containing '{metadata['_table_name']}')."
#         )
#         list_files = list()
#         if isinstance(source_files, list) or Path(source_files).is_file():
#             list_files.append(source_files)
#         else:
#             glob_patterns = [g for g in metadata["_table_name"]]
#             for pattern in glob_patterns:
#                 list_files.extend(
#                     [f for f in Path(source_files).rglob(f"{pattern}*") if f.is_file()]
#                 )
#
#
#
#
#         manager = mp.Manager()
#         errored_files = manager.list()
#         converter_func = partial(
#             _convert_station_file,
#             output_path=rep_nc,
#             errored_files=errored_files,
#             mode=mode,
#             variable_code=variable_code,
#             column_names=column_names,
#             column_dtypes=column_dtypes,
#             **metadata,
#         )
#         with mp.Pool(processes=n_workers) as pool:
#             pool.map(converter_func, list_files)
#             pool.close()
#             pool.join()
#
#
