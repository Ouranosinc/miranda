"""Aggregation module."""
from __future__ import annotations

import logging.config

import xarray as xr
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["aggregations_possible", "aggregate"]

# There needs to be a better way (is there something in xclim?)
_resampling_keys = dict()
_resampling_keys["hour"] = "H"
_resampling_keys["day"] = "D"
_resampling_keys["month"] = "M"
_resampling_keys["year"] = "A"


def aggregations_possible(ds: xr.Dataset, freq: str = "day") -> dict[str, set[str]]:
    """Determine which aggregations are possible based on variables within a dataset.

    Parameters
    ----------
    ds : xarray.Dataset
    freq : str

    Returns
    -------
    dict[str, set[str]]
    """
    logging.info("Determining potential upscaled climate variables.")

    offset, meaning = get_time_frequency(ds, minimum_continuous_period="1h")

    aggregation_legend = dict()
    for v in ["tas", "tdps"]:
        if freq == meaning:
            if not hasattr(ds, v) and (
                hasattr(ds, f"{v}max") and hasattr(ds, f"{v}min")
            ):
                aggregation_legend[f"_{v}"] = {"mean"}
    for variable in ds.data_vars:
        if variable in ["tas", "tdps"]:
            aggregation_legend[variable] = {"max", "mean", "min"}
        elif variable in [
            "evspsblpot",
            "hfls",
            "hfss",
            "hur",
            "hus",
            "pr",
            "prsn",
            "ps",
            "psl",
            "rsds",
            "rss",
            "rlds",
            "rls",
            "snd",
            "snr",
            "snw",
            "swe",
        ]:
            aggregation_legend[variable] = {"mean"}

    return aggregation_legend


def aggregate(ds: xr.Dataset, freq: str = "day") -> dict[str, xr.Dataset]:
    """

    Parameters
    ----------
    ds : xarray.Dataset
    freq : str

    Returns
    -------
    dict[str, xarray.Dataset]
    """
    mappings = aggregations_possible(ds)

    try:
        xarray_agg = _resampling_keys[freq]
    except KeyError:
        xarray_agg = freq

    aggregated = dict()
    for variable, transformations in mappings.items():
        for op in transformations:
            ds_out = xr.Dataset()
            ds_out.attrs = ds.attrs.copy()
            ds_out.attrs["frequency"] = freq

            with xr.set_options(keep_attrs=True):
                if variable.startswith("_"):
                    if op == "mean":
                        var = variable.strip("_")
                        min_var = "".join([var, "min"])
                        max_var = "".join([var, "max"])

                        mean_variable = tas(
                            tasmin=ds[min_var], tasmax=ds[max_var]
                        ).resample(time=xarray_agg)
                        ds_out[var] = mean_variable.mean(dim="time", keep_attrs=True)
                        method = f"time: mean (interval: 1 {freq})"
                        ds_out[var].attrs["cell_methods"] = method
                        aggregated[var] = ds_out
                    continue

                else:
                    if op in {"max", "min"}:
                        transformed = f"{variable}{op}"
                    else:
                        transformed = variable

                    r = ds[variable].resample(time=xarray_agg)
                    ds_out[transformed] = getattr(r, op)(dim="time", keep_attrs=True)
                    method = f"time: {op}{'imum' if op != 'mean' else ''} (interval: 1 {freq})"
                    ds_out[transformed].attrs["cell_methods"] = method
                    aggregated[transformed] = ds_out

    return aggregated
