"""Treatments module."""

from __future__ import annotations
import datetime
import logging

import xarray

from miranda import __version__ as __miranda_version__
from miranda.treatments._dimensions import *
from miranda.treatments._preprocessing import *
from miranda.treatments._variables import *
from miranda.treatments.utils import *
from miranda.units import check_time_frequency


logger = logging.getLogger("miranda.treatments")

VERSION = datetime.datetime.now().strftime("%Y.%m.%d")


def metadata_conversion(d: xarray.Dataset, p: str, m: dict) -> xarray.Dataset:
    """
    Update xarray dataset and data_vars with project-specific metadata fields.

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
    logger.info("Converting metadata to CF-like conventions.")

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
            msg = f"`_miranda_version` not set for project `{p}`. Not appending."
            logger.warning(msg)
    if "_miranda_version" in header:
        del header["_miranda_version"]

    frequency = m["Header"].get("_frequency")
    if frequency:
        if isinstance(frequency, bool):
            _, m["Header"]["frequency"] = check_time_frequency(d)
        elif isinstance(frequency, dict):
            if p in frequency.keys():
                m["Header"]["frequency"] = check_time_frequency(d)
        else:
            logger.warning("`frequency` not set for project. Not appending.")
    if "_frequency" in m["Header"]:
        del m["Header"]["_frequency"]

    # Conditional handling of global attributes based on project name
    for field in [f for f in header.keys() if f.startswith("_")]:
        if isinstance(header[field], list):
            if p in header[field]:
                attr_treatments = header[field][p]
            else:
                msg = f"Attribute handling (`{field}`) not set for project `{p}`. Continuing..."
                logger.warning(msg)
                continue
        elif isinstance(header[field], dict):
            attr_treatments = header[field]
        else:
            msg = f"Attribute treatment configuration for field `{field}` is not properly configured. Verify JSON."
            raise AttributeError(msg)

        if field[1:] in d.attrs:
            msg = f"Overwriting `{field[1:]}` based on JSON configuration."
            logger.warning(msg)
        if field == "_map_attrs":
            for attribute, mapping in attr_treatments.items():
                header[mapping] = d.attrs[attribute]
                del d.attrs[attribute]
        elif field == "_remove_attrs":
            if isinstance(attr_treatments, dict):
                attrs_to_remove = attr_treatments.get(p, [])
            else:
                attrs_to_remove = attr_treatments
            for ff in attrs_to_remove:
                d.attrs.pop(ff, None)
        elif field.startswith("_") and p in attr_treatments:
            header[field[1:]] = attr_treatments[p]
        else:
            header[field[1:]] = attr_treatments
        del header[field]

    # Add global attributes
    d.attrs.update(header)
    d.attrs.update(dict(project=p))

    # Date-based versioning
    if not d.attrs.get("version"):
        d.attrs.update(dict(version=f"v{VERSION}"))

    prev_history = d.attrs.get("history", "")
    history = f"[{datetime.datetime.now()}] Converted variables and modified metadata for CF-like compliance: {prev_history}".strip()
    d.attrs.update(dict(history=history))

    return d
