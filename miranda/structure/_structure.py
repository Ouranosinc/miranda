import logging.config
import os
import shutil
import sys
from pathlib import Path
from types import GeneratorType
from typing import Iterable, List, Mapping, Optional, Union

from miranda.decode import Decoder, guess_project
from miranda.scripting import LOGGING_CONFIG
from miranda.utils import filefolder_iterator
from miranda.validators import GRIDDED_SCHEMA, SIMULATION_SCHEMA, STATION_OBS_SCHEMA

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "build_path_from_schema",
    "structure_datasets",
]


def build_path_from_schema(
    facets: dict, output_folder: Union[str, os.PathLike]
) -> Path:
    """Build a filepath based on a valid data schema.

    Parameters
    ----------
    facets: dict
      Facets for a given dataset.
    output_folder
      Parent folder on which to extend the filetree structure.

    Returns
    -------
    Path
    """
    if facets["type"] == "station-obs":
        STATION_OBS_SCHEMA.validate(facets)
        folder_tree = (
            Path(output_folder)
            / facets["type"]
            / facets["project"]
            / facets["institution"]
            / facets["version"]  # This suggests "date_created"
            / facets["frequency"]
            / facets["variable"]
        )
        if hasattr(facets, "member"):
            return folder_tree / facets["member"]
        return folder_tree

    if facets["type"] in ["forecast", "gridded-obs", "reanalysis"]:
        GRIDDED_SCHEMA.validate(facets)
        return (
            Path(output_folder)
            / facets["type"]
            / facets["institution"]
            / facets["source"]
            / facets["project"]
            / facets["domain"]
            / facets["frequency"]
            / facets["variable"]
        )

    if facets["type"] == "simulation":
        SIMULATION_SCHEMA.validate(facets)
        if facets["project"] == "CORDEX":
            model = facets["driving_model"]
        else:
            model = facets["member"]
        return (
            Path(output_folder)
            / facets["type"]
            / facets["processing_level"]
            / facets["project"]
            / facets["domain"]
            / facets["institution"]
            / facets["source"]
            / model
            / facets["experiment"]
            / facets["member"]
            / facets["frequency"]
            / facets["variable"]
        )

    raise ValueError("No appropriate data schemas found.")


def structure_datasets(
    input_files: Union[str, os.PathLike, List[Union[str, os.PathLike]], GeneratorType],
    output_folder: Union[str, os.PathLike],
    *,
    project: Optional[str] = None,
    guess: bool = True,
    dry_run: bool = False,
    method: str = "copy",
    make_dirs: bool = False,
    filename_pattern: str = "*.nc",
) -> Mapping[Path, Path]:
    """

    Parameters
    ----------
    input_files: str or Path or list of str or Path or GeneratorType
    output_folder: str or Path
    project: {"cordex", "cmip5", "cmip6", "isimip-ft", "reanalysis"}, optional
    guess: bool
      If project not supplied, suggest to decoder that project is the same for all input_files. Default: True.
    dry_run: bool
      Prints changes that would have been made without performing them. Default: False.
    method: {"move", "copy"}
      Method to transfer files to intended location. Default: "move".
    make_dirs:
      Make folder tree if it does not already exist. Default: False.
    filename_pattern: str
      If pattern ends with "zarr", will 'glob' with provided pattern.
      Otherwise, will perform an 'rglob' (recursive) operation.

    Returns
    -------
    dict
    """
    input_files = filefolder_iterator(input_files, filename_pattern)

    if not project and guess:
        for f in input_files:
            project = guess_project(f)
            decoder = Decoder(project)
            decoder.decode(f)
            break
        else:
            raise FileNotFoundError()
        decoder.decode(input_files)
    else:
        decoder = Decoder(project)
        decoder.decode(input_files)

    all_file_paths = dict()
    for file, facets in decoder.file_facets().items():
        output_filepath = build_path_from_schema(facets, output_folder)
        all_file_paths.update({Path(file): output_filepath})

    if make_dirs:
        for new_paths in set(all_file_paths.values()):
            Path(new_paths).mkdir(exist_ok=True, parents=True)

    for file, output_filepath in all_file_paths.items():
        method_mod = ""
        if file.is_dir():
            method_mod = "tree"
        if method.lower() in ["move", "copy"]:
            meth = "Moved" if method.lower() == "move" else "Copied"

            output_file = output_filepath.joinpath(file.name)
            try:
                if not dry_run:
                    if sys.version_info < (3, 9):
                        getattr(shutil, f"{method}{method_mod}")(
                            str(file), str(output_file)
                        )
                    else:
                        getattr(shutil, f"{method}{method_mod}")(file, output_file)
                logging.info(f"{meth} {file.name} to {output_file}.")
            except FileExistsError:
                logging.warning(
                    f"{file.name} already exists at location. Continuing..."
                )

    return all_file_paths
