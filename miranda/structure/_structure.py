import logging
import logging.config
import os

# import shutil
from pathlib import Path
from typing import List, Mapping, Union

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


def _build_path_from_schema(
    schema: dict, output_folder: Union[str, os.PathLike]
) -> Path:
    """Build a filepath based on a validated data schema.

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
            / schema["institute"]
            / schema["project"]
            / schema["version"]  # This suggests "date_created"
            / schema["frequency"]
        )
        if hasattr(schema, "member"):
            return folder_tree / schema["member"]
        else:
            return folder_tree

    elif schema["type"] in ["forecast", "gridded-obs", "reanalysis"]:
        return (
            Path(output_folder)
            / schema["type"]
            / schema["institute"]
            / schema["source"]
            / schema["domain"]
            / schema["frequency"]
        )
    elif schema["type"] == "simulation":
        # FIXME: This currently only works for CORDEX-like data
        return (
            Path(output_folder)
            / schema["type"]
            / schema["processing_level"]
            / schema["project"]
            / schema["domain"]
            / schema["institute"]
            / schema["source"]
            / schema["driving_model"]
            / schema["experiment"]
            / schema["member"]
            / schema["frequency"]
        )


def structure_datasets(
    input_files: Union[str, os.PathLike, List[Union[str, os.PathLike]]],
    output_folder: Union[str, os.PathLike],
    move: bool = False,
) -> Mapping[str, Path]:
    """

    Parameters
    ----------
    input_files: str or Path or list of str or Path
    output_folder: str or Path
    move: bool

    Returns
    -------
    dict
    """
    pass
