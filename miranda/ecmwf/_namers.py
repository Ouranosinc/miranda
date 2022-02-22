import logging
import os
import re
import shutil
from pathlib import Path
from typing import Union

import xarray as xr

__all__ = ["rename_ecmwf_files"]

ECMWF_PROJECT_NAMES = [
    "era5-land",
    "era5-land-monthly-means",
    "era5-pressure-levels",
    "era5-pressure-levels-preliminary-back-extension",
    "era5-single-levels",
    "era5-single-levels-preliminary-back-extension",
]


def rename_ecmwf_files(path: Union[os.PathLike, str]) -> None:
    files = [f for f in Path(path).glob("*.nc")]
    for f in files:
        file_name = str(f.stem)

        ds = xr.open_dataset(f, cache=False)
        var = [d for d in ds.data_vars]
        var_name = str(var[0])

        try:
            x = re.search(r"\d{6}", file_name)
            date = x.group()
        except AttributeError:
            year = int(ds.isel(time=0).time.dt.year)
            month = int(ds.isel(time=0).time.dt.month)
            date = f"{year}{str(month).zfill(2)}"

        names = file_name.split("_")
        projects = [name for name in names if name in ECMWF_PROJECT_NAMES]
        if len(projects) == 1:
            project = projects[0]
        elif len(projects) > 1:
            logging.warning(
                f"More than one project identified for file {f.name}. Verify file naming."
            )
            continue
        else:
            continue

        product = "reanalysis"
        freq = "1hr"
        domain = "NAM"
        institute = "ecmwf"

        new_name_parts = [
            var_name,
            freq,
            institute,
            project,
            product,
            domain,
            date,
        ]
        new_name = f"{'_'.join(new_name_parts)}.nc"
        logging.info(f"Moving {f.name} to {new_name}")

        shutil.move(f, Path(path).joinpath(new_name))
