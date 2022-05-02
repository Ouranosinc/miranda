import functools
import logging
import logging.config
import multiprocessing
import os
import re
import shutil
from datetime import date
from datetime import datetime as dt
from pathlib import Path
from typing import List, Mapping, Optional, Tuple, Union

import xarray as xr

from miranda.gis.subset import subsetting_domains
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["request_era5", "rename_era5_files", "ERA5_PROJECT_NAMES"]


ERA5_PROJECT_NAMES = [
    "era5-land",
    "era5-land-monthly-means",
    "era5-pressure-levels",
    "era5-pressure-levels-preliminary-back-extension",
    "era5-single-levels",
    "era5-single-levels-preliminary-back-extension",
]


def request_era5(
    variables: Optional[Mapping[str, str]],
    projects: List[str],
    *,
    domain: str = "AMNO",
    output_folder: Optional[Union[str, os.PathLike]] = None,
    year_start: Union[str, int] = 1950,
    year_end: Optional[Union[str, int]] = None,
    processes: int = 10,
) -> None:
    """Request ERA5/ERA5-Land from Copernicus Data Store in NetCDF4 format.

    Parameters
    ----------
    variables: Mapping[str, str]
    projects : List[{"era5", "era5-land", "era5-single-levels"}]
    domain : {"GLOBAL", "AMNO", "NAM", "CAN", "QC", "MTL"}
    output_folder : str or os.PathLike, optional
    year_start : int
    year_end : int, optional
    processes : int

    Returns
    -------
    None
    """
    # Variables of interest
    variable_reference = dict()
    variable_reference["era5-land"] = dict(
        tp="total_precipitation",
        v10="10m_v_component_of_wind",
        u10="10m_u_component_of_wind",
        d2m="2m_dewpoint_temperature",
        t2m="2m_temperature",
        pev="potential_evaporation",
        rsn="snow_density",
        sde="snow_depth",
        sd="snow_depth_water_equivalent",
        sf="snowfall",
        swlv1="volumetric_soil_water_layer_1",
        swlv2="volumetric_soil_water_layer_2",
        swlv3="volumetric_soil_water_layer_3",
        swlv4="volumetric_soil_water_layer_4",
    )
    variable_reference[
        "era5", "era-single-levels", "era5-single-levels-preliminary-back-extension"
    ] = dict(
        tp="total_precipitation",
        v10="10m_v_component_of_wind",
        u10="10m_u_component_of_wind",
        d2m="2m_dewpoint_temperature",
        t2m="2m_temperature",
        pev="potential evaporation",
        # sde= Not available for era5
        rsn="snow_density",
        sd="snow_depth",  # note difference in name vs era5-land cf_variable == snw
        sf="snowfall",
        swlv1="volumetric_soil_water_layer_1",
        swlv2="volumetric_soil_water_layer_2",
        swlv3="volumetric_soil_water_layer_3",
        swlv4="volumetric_soil_water_layer_4",
    )

    if year_end is None:
        year_end = date.today().year
    years = range(int(year_start), int(year_end) + 1)

    months = [str(d).zfill(2) for d in range(1, 13)]
    yearmonth = list()
    for y in years:
        for m in months:
            yearmonth.append((y, m))

    project_names = dict()
    if "era5" in projects or "era5-single-levels" in projects:
        project_names["era5-single-levels"] = "reanalysis-era5-single-levels"
    if "era5-land" in projects:
        project_names["era5-land"] = "reanalysis-era5-land"
    if "era5-single-levels-preliminary-back-extension" in projects:
        project_names[
            "era5-single-levels-preliminary-back-extension"
        ] = "reanalysis-era5-single-levels-preliminary-back-extension"

    if output_folder is None:
        target = Path().cwd().joinpath("downloaded")
    else:
        target = output_folder
    Path(target).mkdir(exist_ok=True)
    os.chdir(target)

    for key, p in project_names.items():
        product = p.split("-")[0]
        v_requested = dict()
        variable_reference = next(
            var_list for k, var_list in variable_reference.items() if p in k
        )
        if variables:
            for v in variables:
                if v in variable_reference[key]:
                    v_requested[v] = variable_reference[key][v]
        else:
            v_requested = variable_reference[key]
        proc = multiprocessing.Pool(processes=processes)
        func = functools.partial(_request_direct_era, v_requested, p, domain, product)

        logging.info([func, dt.now().strftime("%Y-%m-%d %X")])

        proc.map(func, yearmonth)
        proc.close()
        proc.join()


def _request_direct_era(
    variables: Mapping[str, str],
    project: str,
    domain: str,
    product: str,
    yearmonth: Tuple[int, str],
):
    """Launch formatted request."""

    try:
        from cdsapi import Client  # noqa
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            f"{_request_direct_era.__name__} requires additional dependencies. "
            "Please install them with `pip install miranda[full]`."
        )

    year, month = yearmonth
    days = [str(d).zfill(2) for d in range(32)]
    times = [f"{str(t).zfill(2)}:00" for t in range(24)]

    if domain.upper() == "AMNO":
        domain = "NAM"
    region = subsetting_domains(domain)

    c = Client()

    if project in ["reanalysis-era5-single-levels", "reanalysis-era5-land"]:
        timestep = "hourly"
    else:
        raise NotImplementedError(project)

    for var in variables.keys():
        netcdf_name = (
            f"{var}_{timestep}_ecmwf_{'-'.join(project.split('-')[1:])}"
            f"_{product}_{domain.upper()}_{year}{month}.nc"
        )

        if Path(netcdf_name).exists():
            logging.info("Dataset %s already exists. Continuing..." % netcdf_name)
            continue

        request_kwargs = dict(
            variable=variables[var],
            year=year,
            month=month,
            day=days,
            time=times,
            area=region,
            format="netcdf",
        )

        if project == "reanalysis-era5-single-levels":
            request_kwargs.update(dict(product_type=product))

        c.retrieve(
            project,
            request_kwargs,
            netcdf_name,
        )


def rename_era5_files(path: Union[os.PathLike, str]) -> None:
    files = [f for f in Path(path).glob("*.nc")]
    for f in files:
        file_name = str(f.stem)

        ds = xr.open_dataset(f, cache=False)
        var = [d for d in ds.data_vars]
        var_name = str(var[0])

        try:
            x = re.search(r"\d{6}", file_name)
            date_found = x.group()
        except AttributeError:
            year = int(ds.isel(time=0).time.dt.year)
            month = int(ds.isel(time=0).time.dt.month)
            date_found = f"{year}{str(month).zfill(2)}"

        names = file_name.split("_")
        projects = [name for name in names if name in ERA5_PROJECT_NAMES]
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
            date_found,
        ]
        new_name = f"{'_'.join(new_name_parts)}.nc"
        logging.info(f"Moving {f.name} to {new_name}")

        shutil.move(f, Path(path).joinpath(new_name))
