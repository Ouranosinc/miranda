from __future__ import annotations

import logging
from typing import Any

from miranda import __version__ as __miranda_version__


def basic_metadata_conversion(
    project: str, metadata: dict
) -> dict[str, dict[str, Any]]:
    """Present basic metadata conversion.

    Parameters
    ----------
    project : str
        Dataset project name.
    metadata : dict
        Metadata definition dictionary for project and variable(s).

    Returns
    -------
    xarray.Dataset
    """
    header = metadata["Header"]

    # Static handling of version global attributes
    miranda_version = header.get("_miranda_version")
    if miranda_version:
        if isinstance(miranda_version, bool):
            header["miranda_version"] = __miranda_version__
        elif isinstance(miranda_version, dict):
            if project in miranda_version.keys():
                header["miranda_version"] = __miranda_version__
        else:
            logging.warning(
                f"`_miranda_version` not set for project `{project}`. Not appending."
            )
    if "_miranda_version" in header:
        del header["_miranda_version"]

    return metadata
