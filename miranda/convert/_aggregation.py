import logging.config
from typing import Dict, Set

import xarray as xr
import xclim.core.options
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["aggregations_possible", "aggregate"]

# There needs to be a better way (is there something in xclim?)
_resampling_keys = dict()
_resampling_keys["hour"] = "H"
_resampling_keys["day"] = "D"
_resampling_keys["month"] = "M"
_resampling_keys["year"] = "A"


def aggregations_possible(ds: xr.Dataset) -> Dict[str, Set[str]]:
    logging.info("Determining potential upscaled climate variables.")

    aggregation_legend = dict()
    for v in ["tas", "tdps"]:
        if not hasattr(ds, v) and (hasattr(ds, f"{v}max") and hasattr(ds, f"{v}min")):
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


def aggregate(ds, freq: str = "day") -> Dict[str, xr.Dataset]:
    mappings = aggregations_possible(ds)

    try:
        xarray_agg = _resampling_keys[freq]
    except KeyError:
        xarray_agg = freq

    aggregated = dict()
    for variable, transformations in mappings.items():
        ds_out = xr.Dataset()
        ds_out.attrs = ds.attrs.copy()
        ds_out.attrs["frequency"] = freq

        with xclim.core.options.set_options(keep_attrs=True):
            if variable.startswith("_"):
                if "mean" in transformations:
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
                for op in transformations:
                    r = ds[variable].resample(time=xarray_agg)
                    ds_out[variable] = getattr(r, op)(dim="time", keep_attrs=True)
                    method = f"time: {op}{'imum' if op != 'mean' else ''} (interval: 1 {freq})"
                    ds_out[variable].attrs["cell_methods"] = method

                    aggregated[variable] = ds_out

    return aggregated
