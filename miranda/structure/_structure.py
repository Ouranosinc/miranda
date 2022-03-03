import logging.config
import os
from pathlib import Path
from typing import List, Mapping, Union

from miranda.decode import (
    decode_ahccd_obs,
    decode_cmip5_name,
    decode_cmip5_netcdf,
    decode_cmip6_name,
    decode_cmip6_netcdf,
    decode_cordex_name,
    decode_cordex_netcdf,
    decode_eccc_obs,
    decode_era5,
    decode_generic_reanalysis,
    decode_isimip_ft_name,
    decode_isimip_ft_netcdf,
    decode_melcc_obs,
    decode_primary_variable,
)
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
        # TODO: Verify whether this is how we want to structure this
        if schema["project"] == "CORDEX":
            model = schema["driving_model"]
        else:
            model = schema["member"]

        return (
            Path(output_folder)
            / schema["type"]
            / schema["processing_level"]
            / schema["project"]
            / schema["domain"]
            / schema["institute"]
            / schema["source"]
            / model
            / schema["experiment"]
            / schema["member"]
            / schema["frequency"]
        )


def structure_datasets(
    input_files: Union[str, os.PathLike, List[Union[str, os.PathLike]]],
    output_folder: Union[str, os.PathLike],
    *,
    project: str,
    move: bool,
    filename_pattern: str = "*",
) -> Mapping[str, Path]:
    """

    Parameters
    ----------
    input_files: str or Path or list of str or Path
    output_folder: str or Path
    project: {"simulation", "reanalysis", "forecast", "gridded-obs, "station-obs"}
    move: bool
    filename_pattern: str

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
