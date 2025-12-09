"""Aggregation module."""

from __future__ import annotations
import logging

import xarray as xr
from xclim.indices import tas

from miranda.units import check_time_frequency


logger = logging.getLogger("miranda.convert.aggregation")

__all__ = ["aggregate", "aggregations_possible"]

# There needs to be a better way (is there something in xclim?)
_resampling_keys = dict()
_resampling_keys["hour"] = "H"
_resampling_keys["day"] = "D"
_resampling_keys["month"] = "M"
_resampling_keys["year"] = "A"


def aggregations_possible(ds: xr.Dataset, freq: str = "day") -> dict[str, set[str]]:
    """
    Determine which aggregations are possible based on variables within a dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset.
    freq : str
        TODO: I'm not entirely certain this is even necessary, but is used to determine whether averages are possible.

    Returns
    -------
    dict
        Mapping of variable names to a set of possible operations (e.g., max, mean, min).

    Notes
    -----
    The function checks first for continuous time periods in the dataset and then
    determines which variables are present and which operations can be performed on them.

    If the dataset has variables that can be aggregated, such as temperature, humidity, and wind speed,
    then the following operations are possible:
    - For temperature: max, mean, min
    - For humidity: max, mean, min
    - For wind speed: max, mean

    For fluxes (e.g., precipitation, evaporation), only the mean operation is available.

    If the dataset has variables that are not present but can be derived (e.g., tas from tasmax and tasmin),
    then the following operations are possible:
    - For derived temperature variables: max, mean, min
    """
    logger.info("Determining potential upscaled climate variables.")

    _, meaning = check_time_frequency(ds, minimum_continuous_period="1H")
    aggregation_legend = {}

    # Variables that are not present in the dataset but that can be derived
    for v in ["tas", "tdps", "hurs"]:
        if freq == meaning:
            if not hasattr(ds, v) and (hasattr(ds, f"{v}max") and hasattr(ds, f"{v}min")):
                aggregation_legend[f"_{v}"] = {"max", "mean", "min"}
                aggregation_legend[f"{v}max"] = {"max", "mean", "min"}
                aggregation_legend[f"{v}min"] = {"max", "mean", "min"}

    # Operations available for variables that are present in the dataset
    for variable in ds.data_vars:
        if variable in ["tas", "ta", "tdps", "tdp", "hurs", "hur", "ts"]:
            aggregation_legend[variable] = {"max", "mean", "min"}
        elif variable in ["sfcWind"]:
            aggregation_legend[variable] = {"max", "mean"}
        # The following variables are expected as fluxes
        elif variable in [
            "CAPE",
            "cfia",
            "evspsblpot",
            "hfls",
            "hfss",
            "huss",
            "hus",
            "pr",
            "prc",
            "prfr",
            "prmod",
            "prra",
            "prrp",
            "prsn",
            "ps",
            "psl",
            "rlds",
            "rls",
            "rsds",
            "rss",
            "snd",
            "sndLand",
            "snr",
            "snw",
            "snwLand",
            "ua",
            "ua100m",
            "uas",
            "va",
            "va100m",
            "vas",
            "winddir",
            "z",
            "zcrd09944",
            "zcrd09975",
            "zcrd10000",
            "20mWind",
            "20mWinddir",
            "40mWind",
        ]:
            aggregation_legend[variable] = {"mean"}

    return aggregation_legend


def aggregate(ds: xr.Dataset, freq: str = "day") -> dict[str, xr.Dataset]:
    """
    Aggregate a dataset to a specified frequency.

    Parameters
    ----------
    ds : xarray.Dataset
    freq : str

    Returns
    -------
    dict[str, xarray.Dataset]
    """
    mappings = aggregations_possible(ds, freq)

    try:
        xarray_agg = _resampling_keys[freq]
    except KeyError:
        xarray_agg = freq

    _ds = ds.copy(deep=True)
    aggregated = {}

    # Calculate the mean variable from max and min variables
    for variable in mappings.copy():
        if variable.startswith("_"):
            var = variable.strip("_")
            min_var = f"{var}min"
            max_var = f"{var}max"
            with xr.set_options(keep_attrs=True):
                _ds[var] = tas(tasmin=ds[min_var], tasmax=ds[max_var])
            offset, meaning = check_time_frequency(ds, minimum_continuous_period="1h")
            method = f"time: mean (interval: {offset} {meaning})"
            _ds[var].attrs["cell_methods"] = method
            del mappings[variable]

    # Aggregate the dataset
    for variable, transformations in mappings.items():
        for op in transformations:
            ds_out = xr.Dataset()
            ds_out.attrs = _ds.attrs.copy()
            ds_out.attrs["frequency"] = freq

            if op in {"max", "min"}:
                transformed = f"{variable}{op}"
            elif op == "mean":
                transformed = variable
            else:
                msg = f"Unsupported operation: {op} for variable {variable}."
                raise ValueError(msg)

            with xr.set_options(keep_attrs=True):
                r = _ds[variable].resample(time=xarray_agg)
                use_snapshot = False
                snapshot_hour = None
                if op == "mean" and "time" in _ds[variable].dims and _ds[variable].time.size > 1:
                    # check if data is only available at a single hour each day
                    other_dims = [d for d in _ds[variable].dims if d != "time"]
                    valid_per_time = _ds[variable].notnull().any(dim=other_dims)
                    valid_by_hour = valid_per_time.groupby(valid_per_time.time.dt.hour).any("time").compute()
                    hours_with_data = valid_by_hour["hour"].values[valid_by_hour.values]
                    if hours_with_data.size == 1:
                        use_snapshot = True
                        snapshot_hour = hours_with_data[0]
                if use_snapshot:
                    # only keep valid data at the detected hour and save as daily values
                    mask_hour = _ds[variable].time.dt.hour == snapshot_hour
                    daily_data = _ds[variable].where(mask_hour, drop=True)
                    ds_out[transformed] = daily_data.resample(time="D").first(keep_attrs=True)
                    ds_out[transformed].attrs["cell_methods"] = "time: point"
                else:
                    ds_out[transformed] = getattr(r, op)(dim="time", keep_attrs=True)
                    method = f"time: {op}{'imum' if op != 'mean' else ''} (interval: 1 {freq})"
                    ds_out[transformed].attrs["cell_methods"] = method
                aggregated[transformed] = ds_out

    return aggregated
