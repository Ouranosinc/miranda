from __future__ import annotations
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from miranda.convert.utils import date_parser


__all__ = [
    "correct_time_entries",
    "correct_var_names",
    "preprocessing_corrections",
]


def correct_time_entries(
    ds: xr.Dataset,
    split: str = "_",
    location: int = -1,
    field: str = "time",
) -> xr.Dataset:
    """
    Correct time entries in dataset.

    Parameters
    ----------
    ds : xarray.Dataset
    split : str
    location : int
    field : str

    Returns
    -------
    xarray.Dataset
    """
    filename = ds.encoding["source"]
    date = date_parser(Path(filename).stem.split(split)[location])
    vals = np.arange(len(ds[field]))
    days_since = f"days since {date}"
    time = xr.coding.times.decode_cf_datetime(vals, units=days_since, calendar="standard")
    ds = ds.assign_coords({field: time})

    prev_history = ds.attrs.get("history", "")
    history = f"Time index recalculated in preprocessing step ({days_since}). {prev_history}"
    ds.attrs.update(dict(history=history))

    return ds


def correct_var_names(ds: xr.Dataset, split: str = "_", location: int = 0) -> xr.Dataset:
    """
    Correct variable names in dataset.

    Parameters
    ----------
    ds : xarray.Dataset
    split : str
    location : int

    Returns
    -------
    xarray.Dataset
    """
    filename = ds.encoding["source"]
    new_name = Path(filename).stem.split(split)[location]
    old_name = list(ds.data_vars.keys())[0]

    prev_history = ds.attrs.get("history", "")
    history = f"Variable renamed in preprocessing step ({old_name}: {new_name}). {prev_history}"
    ds.attrs.update(dict(history=history))

    return ds.rename({old_name: new_name})


def preprocessing_corrections(ds: xr.Dataset, configuration: dict[str, Any]) -> xr.Dataset:
    """
    Corrections function dispatcher to ensure minimal dataset validity on open.

    Parameters
    ----------
    ds : xarray.Dataset
    configuration : dict

    Returns
    -------
    xarray.Dataset
    """

    def _preprocess_correct(d: xr.Dataset, *, ops: list[partial]) -> xr.Dataset:
        for correction in ops:
            d = correction(d)
        return d

    correction_fields = configuration.get("_preprocess")
    if correction_fields:
        preprocess_ops = []
        for field in correction_fields:
            if field == "_variable_name":
                preprocess_ops.append(partial(correct_var_names, **correction_fields[field]))
            if field == "_time":
                preprocess_ops.append(partial(correct_time_entries, **correction_fields[field]))
        if preprocess_ops:
            corrector = partial(_preprocess_correct, ops=preprocess_ops)
            return corrector(ds)
    return ds
