from __future__ import annotations

import datetime
import functools
import logging.config
import multiprocessing
import os
import re
import shutil
import warnings
from collections.abc import Sequence
from datetime import datetime as dt
from pathlib import Path
from typing import Any

import xarray as xr

try:
    from cdsapi import Client  # noqa
except ModuleNotFoundError:
    warnings.warn(
        "ERA5 downloading capabilities requires additional dependencies. "
        "Please install them with `pip install miranda[remote]`."
    )

from miranda.gis import subsetting_domains
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

logging.config.dictConfig(LOGGING_CONFIG)


__all__ = ["request_era5", "rename_era5_files", "ERA5_PROJECT_NAMES"]


ERA5_PROJECT_NAMES = [
    "era5-land",
    "era5-land-monthly-means",
    "era5-pressure-levels",
    "era5-pressure-levels-monthly-means",
    "era5-pressure-levels-preliminary-back-extension",
    "era5-pressure-levels-monthly-means-preliminary-back-extension",
    "era5-single-levels",
    "era5-single-levels-monthly-means",
    "era5-single-levels-preliminary-back-extension",
    "era5-single-levels-monthly-means-preliminary-back-extension",
]


def request_era5(
    projects: str | list[str],
    *,
    variables: str | Sequence[str] | None = None,
    domain: str = "AMNO",
    pressure_levels: list[int] | None = None,
    separate_pressure_levels: bool = True,
    output_folder: str | os.PathLike | None = None,
    year_start: str | int | None = None,
    year_end: str | int | None = None,
    dry_run: bool = False,
    processes: int = 10,
    url: str | None = None,
    key: str | None = None,
) -> None:
    """Request ERA5/ERA5-Land from Copernicus Data Store in NetCDF4 format.

    Parameters
    ----------
    projects : str or List[str]
        Allowed keys: {"era5-land", "era5-land-monthly-means", "era5-single-levels", "era5-single-levels-monthly-means",
        "era5-single-levels-preliminary-back-extension", "era5-single-levels-monthly-means-preliminary-back-extension",
        "era5-pressure-levels", "era5-pressure-levels-monthly-means", "era5-pressure-levels-preliminary-back-extension",
        "era5-pressure-levels-monthly-means-preliminary-back-extension"}
    variables : str or Sequence[str]
        Variable codes requested. If None, will attempt all hardcoded variables supported by miranda converter.
    domain : {"GLOBAL", "AMNO", "NAM", "CAN", "QC", "MTL"}
        Geographic domain requested. Default: "AMNO" (North America).
    pressure_levels : List[int], optional
        If set and project requested has pressure levels, will download specific pressure levels.
    separate_pressure_levels : bool
        Whether to separate files for each pressure level. Default: True
    output_folder : str or os.PathLike, optional
        Folder to send files to. If None, will create a "downloaded" folder in current working directory.
    year_start : int, optional
        Starting year for data download. If None, will download from first available year for project.
    year_end : int, optional
        End year for data download. If None, will download files for current year and two months prior to present day.
    dry_run : bool
        Do not send request. For debugging purposes.
    processes : int
        The number of simultaneous download requests. Default: 10.
    url : str, optional
        URL for Copernicus Data Store API (if not already using .cdsapirc)
    key : str, optional
      Personal access key for Copernicus Data Store (if not already using .cdsapirc)

    Returns
    -------
    None
    """
    # Variables of interest
    variable_reference = dict()
    variable_reference["era5-land", "era5-land-monthly-means"] = dict(
        d2m="2m_dewpoint_temperature",
        pev="potential_evaporation",
        rsn="snow_density",
        sd="snow_depth_water_equivalent",
        sde="snow_depth",
        sf="snowfall",
        slhf="surface_latent_heat_flux",
        sp="surface_pressure",
        sshf="surface_sensible_heat_flux",
        ssr="surface_net_solar_radiation",
        ssrd="surface_solar_radiation_downwards",
        str="surface_net_thermal_radiation",
        strd="surface_thermal_radiation_downwards",
        swvl1="volumetric_soil_water_layer_1",
        swvl2="volumetric_soil_water_layer_2",
        swvl3="volumetric_soil_water_layer_3",
        swvl4="volumetric_soil_water_layer_4",
        t2m="2m_temperature",
        tp="total_precipitation",
        u10="10m_u_component_of_wind",
        v10="10m_v_component_of_wind",
    )

    era5_single_levels = variable_reference[
        "era5-land", "era5-land-monthly-means"
    ].copy()
    del era5_single_levels["sde"]  # sde is not available for era5
    era5_single_levels.update(
        msl="mean_sea_level_pressure", sd="snow_depth", ptype="precipitation_type"
    )  # note difference in name vs era5-land cf_variable == snw"
    variable_reference[
        "era5-single-levels",
        "era5-single-levels-monthly-means",
        "era5-single-levels-preliminary-back-extension",
        "era5-single-levels-monthly-means-preliminary-back-extension",
    ] = era5_single_levels

    variable_reference[
        "era5-pressure-levels",
        "era5-pressure-levels-monthly-means",
        "era5-pressure-levels-preliminary-back-extension",
        "era5-pressure-levels-monthly-means-preliminary-back-extension",
    ] = dict(
        cc="fraction_of_cloud_cover",
        r="relative_humidity",
        q="specific_humidity",
        t="temperature",
        u="u_component_of_wind",
        v="v_component_of_wind",
        z="geopotential",
    )

    if output_folder is None:
        target = Path().cwd().joinpath("downloaded")
    else:
        target = output_folder
    Path(target).mkdir(exist_ok=True)
    os.chdir(target)

    project_names = dict()
    if isinstance(projects, str):
        projects = [projects]
    for project in projects:
        project_names[project] = f"reanalysis-{project}"

    for project_name, request_code in project_names.items():
        if year_start is None:
            if "preliminary-back-extension" in project_name or project_name.startswith(
                "era5-land"
            ):
                project_year_start = 1940
            else:
                project_year_start = 1959
        else:
            project_year_start = year_start

        if year_end is None:
            if "preliminary-back-extension" in project_name:
                project_year_end = 1978
            else:
                project_year_end = dt.today().year
        else:
            project_year_end = year_end

        years = range(int(project_year_start), int(project_year_end) + 1)

        # Allow at least a two-month lag from current month before attempting to collect data
        months = [str(d).zfill(2) for d in range(1, 13)]
        yearmonth = list()

        today = datetime.date.today()
        two_months_ago = today - datetime.timedelta(60)
        for y in years:
            if "monthly-means" in project_name:
                request_date = datetime.date(y, today.month, today.day)
                if y < today.year:
                    yearmonth.append((y, months))
                elif two_months_ago.year == today.year:
                    yearmonth.append(
                        (y, [str(d).zfill(2) for d in range(1, two_months_ago.month)])
                    )
            else:
                for m in months:
                    request_date = datetime.date(y, int(m), today.day)
                    if y < today.year - 2:
                        yearmonth.append((y, m))
                    else:
                        if request_date < two_months_ago:
                            yearmonth.append((y, m))

        if "monthly-means" in project_name:
            product = "monthly_averaged_reanalysis"
        else:
            product = "reanalysis"

        v_requested = dict()
        try:
            variable_reference = next(
                var_list
                for k, var_list in variable_reference.items()
                if project_name in k
            )
        except StopIteration:
            return

        if variables:
            if isinstance(variables, str):
                variables = [variables]
            for v in variables:
                if v in variable_reference:
                    v_requested[v] = variable_reference[v]
        else:
            v_requested = variable_reference

        if "pressure-levels" in project_name:
            pressure_levels_requested = [str(i) for i in pressure_levels]
        else:
            pressure_levels_requested = None

        if not dry_run:
            client_kwargs = dict()
            if url:
                client_kwargs["url"] = url
            if key:
                client_kwargs["key"] = key

            client = Client(**client_kwargs)
        else:
            client = None

        proc = multiprocessing.Pool(processes=processes)
        func = functools.partial(
            _request_direct_era,
            v_requested,
            request_code,
            domain,
            pressure_levels_requested,
            separate_pressure_levels,
            product,
            dry_run,
            client,
        )

        logging.info([func, dt.now().strftime("%Y-%m-%d %X")])

        proc.map(func, yearmonth)
        proc.close()
        proc.join()


def _request_direct_era(
    variables: dict[str, str],
    project: str,
    domain: str,
    pressure_levels: list[str] | None,
    separate_pressure_level_requests: bool,
    product: str,
    dry_run: bool,
    client: Any | None,
    yearmonth: tuple[int, str],
):
    """Launch formatted request."""

    def __request(nc_name: str, p: str, rq_kwargs: dict[str, str], c: Any):
        if Path(nc_name).exists():
            logging.info(f"Dataset {nc_name} already exists. Continuing...")
            return

        if not dry_run:
            c.retrieve(
                p,
                rq_kwargs,
                nc_name,
            )
        else:
            logging.info(p)
            logging.info(rq_kwargs)
            logging.info(nc_name)

    year, month = yearmonth
    year_month = f"{year}{month if not isinstance(month, list) else ''}"

    if "monthly-means" in project:
        times = "00:00"
        day = dict()
        timestep = "mon"
    else:
        times = [f"{str(t).zfill(2)}:00" for t in range(24)]
        day = dict(day=[str(d).zfill(2) for d in range(32)])
        timestep = "1h"

    if domain.upper() == "AMNO":
        domain = "NAM"

    region = subsetting_domains(domain)

    for var in variables.keys():
        request_kwargs = dict(
            variable=variables[var],
            year=year,
            month=month,
            **day,
            time=times,
            area=region,
            format="netcdf",
        )

        if (
            "reanalysis-era5-single-levels" in project
            or "reanalysis-era5-pressure-levels" in project
            or "monthly-means" in project
        ):
            request_kwargs.update(dict(product_type=product))

        if pressure_levels:
            if separate_pressure_level_requests:
                for level in pressure_levels:
                    request_kwargs.update(dict(pressure_level=[level]))
                    netcdf_name = (
                        f"{var}{level}_{timestep}_ecmwf_{'-'.join(project.split('-')[1:])}"
                        f"_{product.replace('_', '-')}_{domain.upper()}_{year_month}.nc"
                    )
                    __request(netcdf_name, project, request_kwargs)
                continue
            else:
                request_kwargs.update(dict(pressure_level=pressure_levels))

        netcdf_name = (
            f"{var}_{timestep}_ecmwf_{'-'.join(project.split('-')[1:])}"
            f"_{product.replace('_', '-')}_{domain.upper()}_{year_month}.nc"
        )
        __request(netcdf_name, project, request_kwargs, client)


def rename_era5_files(path: os.PathLike | str) -> None:
    """Rename badly named ERA5 files.

    Notes
    -----
    Requires that the proper ERA5 project name is in the filename, separated by underscores.
    Assumes that the data

    Parameters
    ----------
    path: os.PathLike or str
      Path to a folder containing netcdf files

    Returns
    -------
    None

    """
    files = Path(path).glob("*.nc")
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

        try:
            freq_parts = get_time_frequency(ds)
            freq = f"{freq_parts[0]}{freq_parts[1]}"
        except ValueError:
            logging.error(
                f"Unable to parse the time frequency for variable `{var_name}` "
                f"in file `{f.name}`. Verify data integrity before retrying."
            )
            continue

        names = file_name.split("_")
        projects = [name for name in names if name in ERA5_PROJECT_NAMES]
        if len(projects) == 1:
            project = projects.pop()
        elif len(projects) > 1:
            logging.warning(
                f"More than one project identified for file {f.name}. Verify file naming."
            )
            continue
        else:
            logging.warning("No project string found in filename.")
            continue

        product = "reanalysis"
        institute = "ecmwf"

        new_name_parts = [
            var_name,
            freq,
            institute,
            project,
            product,
            date_found,
        ]
        new_name = f"{'_'.join(new_name_parts)}.nc"
        logging.info(f"Moving {f.name} to {new_name}")

        shutil.move(f, Path(path).joinpath(new_name))
