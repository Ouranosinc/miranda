from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_config_folder = Path(__file__).resolve().parent / "configs"


__all__ = ["load_json_data_mappings"]


def load_json_data_mappings(project: str) -> dict[str, Any]:
    """Load JSON mappings for supported dataset conversions.

    Parameters
    ----------
    project : str

    Returns
    -------
    dict[str, Any]
    """
    if project == "eccc-homogenized":
        metadata_definition = json.load(
            open(_config_folder / "eccc-homogenized_attrs.json")
        )
    elif project == "eccc-obs":
        metadata_definition = json.load(open(_config_folder / "eccc-obs_attrs.json"))
    elif project == "eccc-obs-summary":
        metadata_definition = json.load(
            open(_config_folder / "eccc-obs-summary_attrs.json")
        )
    else:
        raise NotImplementedError(f"Project not supported: {project}")

    return metadata_definition
