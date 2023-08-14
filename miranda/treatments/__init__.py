"""Treatments module."""

from __future__ import annotations

import datetime
import logging.config

import xarray as xr

from miranda import __version__ as __miranda_version__
from miranda.scripting import LOGGING_CONFIG
from miranda.treatments._dimensions import *
from miranda.treatments._preprocessing import *
from miranda.treatments._variables import *
from miranda.treatments.utils import *
from miranda.units import get_time_frequency

logging.config.dictConfig(LOGGING_CONFIG)
VERSION = datetime.datetime.now().strftime("%Y.%m.%d")


def metadata_conversion(d: xr.Dataset, p: str, m: dict) -> xr.Dataset:
    """Update xarray dataset and data_vars with project-specific metadata fields.

    Parameters
    ----------
    d : xarray.Dataset
        Dataset with metadata to be updated.
    p : str
        Dataset project name.
    m : dict
        Metadata definition dictionary for project and variable(s).

    Returns
    -------
    xarray.Dataset
    """
    logging.info("Converting metadata to CF-like conventions.")

    header = m["Header"]

    # Static handling of version global attributes
    miranda_version = header.get("_miranda_version")
    if miranda_version:
        if isinstance(miranda_version, bool):
            header["miranda_version"] = __miranda_version__
        elif isinstance(miranda_version, dict):
            if p in miranda_version.keys():
                header["miranda_version"] = __miranda_version__
        else:
            logging.warning(
                f"`_miranda_version` not set for project `{p}`. Not appending."
            )
    if "_miranda_version" in header:
        del header["_miranda_version"]

    frequency = m["Header"].get("_frequency")
    if frequency:
        if isinstance(frequency, bool):
            _, m["Header"]["frequency"] = get_time_frequency(d)
        elif isinstance(frequency, dict):
            if p in frequency.keys():
                m["Header"]["frequency"] = get_time_frequency(d)
        else:
            logging.warning("`frequency` not set for project. Not appending.")
    if "_frequency" in m["Header"]:
        del m["Header"]["_frequency"]

    # Conditional handling of global attributes based on project name
    for field in [f for f in header if f.startswith("_")]:
        if isinstance(header[field], list):
            if p in header[field]:
                attr_treatment = header[field][p]
            else:
                logging.warning(
                    f"Attribute handling (`{field}`) not set for project `{p}`. Continuing..."
                )
                continue
        elif isinstance(header[field], dict):
            attr_treatment = header[field]
        else:
            raise AttributeError(
                f"Attribute treatment configuration for field `{field}` is not properly configured. Verify JSON."
            )

        if field == "_map_attrs":
            for attribute, mapping in attr_treatment.items():
                header[mapping] = d.attrs[attribute]
                del d.attrs[attribute]
        elif field == "_remove_attrs":
            for ff in attr_treatment:
                del d.attrs[ff]
        else:
            if field[1:] in d.attrs:
                logging.warning(
                    f"Overwriting `{field[1:]}` based on JSON configuration."
                )
            header[field[1:]] = attr_treatment

        del header[field]

    # Add global attributes
    d.attrs.update(header)
    d.attrs.update(dict(project=p))

    # Date-based versioning
    if not d.attrs.get("version"):
        d.attrs.update(dict(version=f"v{VERSION}"))

    prev_history = d.attrs.get("history", "")
    history = (
        f"[{datetime.datetime.now()}] "
        "Converted variables and modified metadata for CF-like compliance: "
        f"{prev_history}".strip()
    )
    d.attrs.update(dict(history=history))

    return d
