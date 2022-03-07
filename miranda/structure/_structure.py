import logging.config
import os
import shutil
from pathlib import Path
from typing import List, Mapping, Optional, Union

import schema

from miranda.decode import Decoder, guess_project
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

STATION_OBS_SCHEMA = schema.Schema(
    {
        "project": str,
        "institution": str,
        "source": str,
        "frequency": str,
        "variable": str,
        "version": str,
        "type": "station-obs",
    },
    ignore_extra_keys=True,
)


GRIDDED_SCHEMA = schema.Schema(
    {
        "project": str,
        "institution": str,
        "source": str,
        "domain": str,
        "frequency": str,
        "variable": str,
        "type": schema.And(
            str,
            lambda f: f in ["forecast", "gridded-obs", "reanalysis"],
        ),
    },
    ignore_extra_keys=True,
)

SIMULATION_SCHEMA = schema.Schema(
    {
        "type": "simulation",
        "processing_level": schema.And(str, lambda f: f in ["raw", "biasadjusted"]),
        "project": str,
        "institution": str,
        "source": str,
        "domain": str,
        schema.Or("driving_model", "member"): str,
        "experiment": str,
        "frequency": str,
        "member": str,
        "variable": str,
    },
    ignore_extra_keys=True,
)


def _build_path_from_schema(
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
        else:
            return folder_tree

    elif facets["type"] in ["forecast", "gridded-obs", "reanalysis"]:
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
    elif facets["type"] == "simulation":
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


def structure_datasets(
    input_files: Union[str, os.PathLike, List[Union[str, os.PathLike]]],
    output_folder: Union[str, os.PathLike],
    *,
    project: Optional[str] = None,
    guess: bool = True,
    dry_run: bool = True,
    copy: bool = False,
    make_dirs: bool = False,
    filename_pattern: str = "*.nc",
) -> Mapping[Path, Path]:
    """

    Parameters
    ----------
    input_files: str or Path or list of str or Path
    output_folder: str or Path
    project: {"cordex", "cmip5", "cmip6", "isimip-ft"}, optional
    guess: bool
      If project not supplied, suggest to decoder that project is the same for all input_files. Default: True.
    dry_run: bool
      Prints changes that would have been made without performing them. Default: True.
    copy: bool
      Make a copy of files to intended location. Default: False.
    make_dirs:
      Make folder tree if it does not already exist. Default: False.
    filename_pattern: str
      If given a path to a directory, will 'glob' with provided pattern.

    Returns
    -------
    dict
    """
    if isinstance(input_files, (Path, str)):
        input_files = Path(input_files)
        if input_files.is_dir():
            input_files = sorted(list(input_files.glob(filename_pattern)))
    elif isinstance(input_files, list):
        input_files = sorted(Path(p) for p in input_files)
    else:
        raise NotImplementedError()

    if not project and guess:
        project = guess_project(input_files[0])

    decoder = Decoder(project)
    decoder.decode(input_files)

    all_file_paths = dict()
    for file, facets in decoder.file_facets().items():
        output_filepath = _build_path_from_schema(facets, output_folder)
        all_file_paths.update({Path(file): output_filepath})

    if make_dirs:
        for new_paths in set(all_file_paths.values()):
            Path(new_paths).mkdir(exist_ok=True, parents=True)

    for file, output_filepath in all_file_paths.items():
        if copy:
            logging.info(f"Copied {file} to {output_filepath}.")
            if not dry_run:
                shutil.copy(file, output_filepath)
        else:
            logging.info(f"Moved {file} to {output_filepath}.")
            if not dry_run:
                shutil.move(file, output_filepath)

    return all_file_paths
