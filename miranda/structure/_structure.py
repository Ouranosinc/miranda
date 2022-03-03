import logging.config
import os
import shutil
from pathlib import Path
from typing import List, Mapping, Optional, Union

from miranda.decode import Decoder, guess_project
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


def _build_path_from_schema(
    schema: dict, output_folder: Union[str, os.PathLike]
) -> Path:
    """Build a filepath based on a valid data schema.

    Parameters
    ----------
    schema: dict
      Validated facet schema for a given dataset.
    output_folder
      Parent folder on which to extend the filetree structure.

    Returns
    -------
    Path
    """
    if schema["type"] == "station-obs":
        folder_tree = (
            Path(output_folder)
            / schema["type"]
            / schema["project"]
            / schema["institution"]
            / schema["project"]
            / schema["version"]  # This suggests "date_created"
            / schema["frequency"]
            / schema["variable"]
        )
        if hasattr(schema, "member"):
            return folder_tree / schema["member"]
        else:
            return folder_tree

    elif schema["type"] in ["forecast", "gridded-obs", "reanalysis"]:
        return (
            Path(output_folder)
            / schema["type"]
            / schema["project"]
            / schema["institution"]
            / schema["source"]
            / schema["domain"]
            / schema["frequency"]
            / schema["variable"]
        )
    elif schema["type"] == "simulation":
        # TODO: Verify whether this is how we want to structure this
        if schema["project"] == "CORDEX":
            model = schema["driving_model"]
        else:
            model = schema["member"]

        return (
            Path(output_folder)
            / schema["type"]
            / schema["project"]
            / schema["processing_level"]
            / schema["project"]
            / schema["domain"]
            / schema["institution"]
            / schema["source"]
            / model
            / schema["experiment"]
            / schema["member"]
            / schema["frequency"]
            / schema["variable"]
        )


def structure_datasets(
    input_files: Union[str, os.PathLike, List[Union[str, os.PathLike]]],
    output_folder: Union[str, os.PathLike],
    *,
    project: Optional[str] = None,
    guess: bool = True,
    copy: bool = False,
    filename_pattern: str = "*.nc",
) -> Mapping[str, Path]:
    """

    Parameters
    ----------
    input_files: str or Path or list of str or Path
    output_folder: str or Path
    project: {"cordex", "cmip5", "cmip6", "isimip-ft"}, optional
    guess: bool
      If project not supplied, suggest to decoder that project is the same for all input_files. Default: True.
    copy: bool
      Make a copy of files to intended location. Default: False.
    filename_pattern: str

    Returns
    -------
    dict
    """
    if isinstance(input_files, (Path, str)):
        input_files = Path(input_files)
        if input_files.is_dir():
            input_files = sorted(list(input_files.rglob(filename_pattern)))
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
        all_file_paths.update({Path(file).name: output_filepath})

        if copy:
            shutil.copy(file, output_filepath)

    return all_file_paths
