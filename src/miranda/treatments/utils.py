"""Utility functions for GIS operations."""

from __future__ import annotations
import inspect
import json
from pathlib import Path
from typing import Any


__all__ = [
    "load_json_data_mappings",
]


def _get_section_entry_key(meta: dict, entry: str, var: str, key: str, project: str) -> Any:
    """
    Get a specific key from a section of the metadata.

    Parameters
    ----------
    meta : dict
        The metadata dictionary.
    entry : str
        The entry to look for.
    var : str
        The variable to look for.
    key : str
        The key to look for.
    project : str
        The project name.

    Returns
    -------
    Any
        The value of the key.
    """
    var_meta = meta[entry].get(var, {})
    if key in var_meta:
        if isinstance(var_meta[key], dict):
            config = var_meta[key].get(project)
            if config is None and "all" in var_meta[key].keys():
                config = var_meta[key].get("all")
            return config
        return var_meta[key]
    return None


def _iter_entry_key(ds, meta, entry, key, project) -> tuple[str, Any]:
    """
    Iterate through entry keys.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    meta : dict
        The metadata dictionary.
    entry : str
        The entry to look for.
    key : str
        The key to look for.
    project : str
        The project name.

    Yields
    ------
    tuple[str, Any]
        The variable and value.
    """
    for vv in set(ds.data_vars).intersection(meta[entry]):
        val = _get_section_entry_key(meta, entry, vv, key, project)
        yield vv, val


def load_json_data_mappings(project: str, configurations: dict[str, Path] | None = None) -> dict[str, Any]:
    """
    Load JSON mappings for supported dataset conversions.

    Parameters
    ----------
    project : str
        The project name.
    configurations : dict, optional
        Configuration files for the project.
        If not provided, the function will try to find the configuration files in the `configs` folder.

    Returns
    -------
    dict[str, Any]
        The metadata definition.
    """
    if configurations is None:
        calling_frame = inspect.currentframe().f_back
        calling_file_path = calling_frame.f_globals["__file__"]
        config_folder = Path(calling_file_path).parent / "configs"

        configurations = {}
        for configuration in config_folder.glob("*attrs.json"):
            project_config = str(configuration.stem).split("_")[0]
            if "|" in project:
                for p in project_config.split("|"):
                    configurations[p] = configuration
            configurations[project_config] = configuration

    if project in configurations.keys():
        config_file = configurations[project]
        metadata_definition = json.load(config_file.open())
        return metadata_definition
    else:
        raise NotImplementedError(f"Project not supported: {project}")
