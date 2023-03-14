import hashlib
import logging.config
import os
import re
from pathlib import Path
from typing import Dict, Union

import xarray as xr
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "daily_aggregation",
    "find_version_hash",
]


def daily_aggregation(ds: xr.Dataset, keys_only: bool = False) -> Dict[str, xr.Dataset]:
    logging.info("Creating daily upscaled climate variables.")

    daily_dataset = dict()
    for variable in ds.data_vars:
        if variable in ["tas", "tdps"]:
            # Some looping to deal with memory consumption issues
            for v, func in {
                f"{variable}max": "max",
                f"{variable}min": "min",
                f"{variable}": "mean",
            }.items():
                if not keys_only:
                    ds_out = xr.Dataset()
                    ds_out.attrs = ds.attrs.copy()
                    ds_out.attrs["frequency"] = "day"

                    method = f"time: {func}{'imum' if func != 'mean' else ''} (interval: 1 day)"
                    ds_out.attrs["cell_methods"] = method

                    if v == "tas" and not hasattr(ds, "tas"):
                        ds_out[v] = tas(tasmax=ds.tasmax, tasmin=ds.tasmin)
                    else:
                        # Thanks for the help, xclim contributors
                        r = ds[variable].resample(time="D")
                        ds_out[v] = getattr(r, func)(dim="time", keep_attrs=True)

                    daily_dataset[v] = ds_out
                    del ds_out
                else:
                    daily_dataset[v] = []

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
            if not keys_only:
                ds_out = xr.Dataset()
                ds_out.attrs = ds.attrs.copy()
                ds_out.attrs["frequency"] = "day"
                ds_out.attrs["cell_methods"] = "time: mean (interval: 1 day)"
                logging.info(f"Converting {variable} to daily time step (daily mean).")
                ds_out[variable] = (
                    ds[variable].resample(time="D").mean(dim="time", keep_attrs=True)
                )

                daily_dataset[variable] = ds_out
                del ds_out
            else:
                daily_dataset[variable] = []
        else:
            continue

    return daily_dataset


def find_version_hash(file: Union[os.PathLike, str]) -> Dict:
    """

    Parameters
    ----------
    file : Path or str

    Returns
    -------
    dict
    """

    def _get_hash(f):
        hash_sha256_writer = hashlib.sha256()
        with open(f, "rb") as f_opened:
            hash_sha256_writer.update(f_opened.read())
        sha256sum = hash_sha256_writer.hexdigest()
        logging.info(f"Calculated sha256sum (starting: {sha256sum[:6]})")
        del hash_sha256_writer
        return sha256sum

    version_info = dict()
    possible_version = Path(file).parent.name
    if re.match(r"^v\d+", possible_version, re.IGNORECASE):
        version_info["version"] = Path(file).parent.name
        version_info["sha256sum"] = _get_hash(file)

    else:
        file_identity = str(Path(file).name).split(".")[0]
        possible_version_signature = Path(file).parent.glob(f"{file_identity}.*")
        for sig in possible_version_signature:
            found_version = re.search(r"\.(v\d+.+)$", sig.name, re.IGNORECASE)
            if found_version:
                try:
                    version_info["version"] = found_version.group()
                    version_info["sha256sum"] = int(sig.open().read())
                except ValueError:
                    continue
                break
        else:
            version_info["version"] = "vNotFound"
            version_info["sha256sum"] = _get_hash(file)

    return version_info
