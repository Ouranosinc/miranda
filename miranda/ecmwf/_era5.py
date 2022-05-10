import datetime
import functools
import logging
import logging.config
import multiprocessing
import os
import re
import shutil
from datetime import datetime as dt
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple, Union

import xarray as xr

from miranda.gis.subset import subsetting_domains
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

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
    projects: List[str],
    *,
    variables: Optional[Mapping[str, str]] = None,
    domain: str = "AMNO",
    pressure_levels: Optional[List[int]] = None,
    separate_pressure_levels: bool = True,
    output_folder: Optional[Union[str, os.PathLike]] = None,
    year_start: Optional[Union[str, int]] = None,
    year_end: Optional[Union[str, int]] = None,
    dry_run: bool = False,
    processes: int = 10,
) -> None:
    """Request ERA5/ERA5-Land from Copernicus Data Store in NetCDF4 format.

    Parameters
    ----------
    variables: Mapping[str, str]
    projects : List[{"era5-land", "era5-single-levels", "era5-single-levels-preliminary-back-extension", "era5-pressure-levels",  "era5-pressure-levels-preliminary-back-extension"}]
    domain : {"GLOBAL", "AMNO", "NAM", "CAN", "QC", "MTL"}
    pressure_levels: List[int], optional
    separate_pressure_levels: bool
      Separate files for each pressure level. Default: True
    output_folder : str or os.PathLike, optional
    year_start : int, optional
    year_end : int, optional
    dry_run: bool
      Do not send request. For debugging purposes.
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
        sp="surface_pressure",
        sshf="surface_sensible_heat_flux",
        slhf="surface_latent_heat_flux",
        ssrd="surface_solar_radiation_downwards",
        strd="surface_thermal_radiation_downwards",
        swlv1="volumetric_soil_water_layer_1",
        swlv2="volumetric_soil_water_layer_2",
        swlv3="volumetric_soil_water_layer_3",
        swlv4="volumetric_soil_water_layer_4",
    )
    variable_reference[
        "era5-single-levels", "era5-single-levels-preliminary-back-extension"
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
        sp="surface_pressure",
        sshf="surface_sensible_heat_flux",
        slhf="surface_latent_heat_flux",
        ssrd="surface_solar_radiation_downwards",
        strd="surface_thermal_radiation_downwards",
        swlv1="volumetric_soil_water_layer_1",
        swlv2="volumetric_soil_water_layer_2",
        swlv3="volumetric_soil_water_layer_3",
        swlv4="volumetric_soil_water_layer_4",
    )
    variable_reference[
        "era5-pressure-levels", "era5-pressure-levels-preliminary-back-extension"
    ] = dict(z="geopotential")

    if output_folder is None:
        target = Path().cwd().joinpath("downloaded")
    else:
        target = output_folder
    Path(target).mkdir(exist_ok=True)
    os.chdir(target)

    project_names = dict()
    for project in projects:
        project_names[project] = f"reanalysis-{project}"

    for project_name, request_code in project_names.items():
        if year_start is None:
            if "back-extension" in project_name or project_name == "era5-land":
                project_year_start = 1950
            else:
                project_year_start = 1979
        else:
            project_year_start = year_start

        if year_end is None:
            if "back-extension" in project_name:
                project_year_end = 1978
            else:
                project_year_end = dt.today().year
        else:
            project_year_end = year_end

        years = range(int(project_year_start), int(project_year_end) + 1)

        months = [str(d).zfill(2) for d in range(1, 13)]
        yearmonth = list()
        for y in years:
            for m in months:
                request_date = datetime.date(y, int(m), 1)
                two_months_ago = datetime.date.today() - datetime.timedelta(60)
                if request_date < two_months_ago:
                    yearmonth.append((y, m))

        product = request_code.split("-")[0]
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
            for v in variables:
                if v in variable_reference:
                    v_requested[v] = variable_reference[v]
        else:
            v_requested = variable_reference

        if "pressure-levels" in project_name:
            pressure_levels_requested = [str(i) for i in pressure_levels]
        else:
            pressure_levels_requested = None

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
        )

        logging.info([func, dt.now().strftime("%Y-%m-%d %X")])

        proc.map(func, yearmonth)
        proc.close()
        proc.join()


def _request_direct_era(
    variables: Mapping[str, str],
    project: str,
    domain: str,
    pressure_levels: Optional[List[str]],
    separate_pressure_level_requests: bool,
    product: str,
    dry_run: bool,
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

    if "monthly-means" in project:
        raise NotImplementedError(project)
    timestep = "1h"

    for var in variables.keys():
        request_kwargs = dict(
            variable=variables[var],
            year=year,
            month=month,
            day=days,
            time=times,
            area=region,
            format="netcdf",
        )

        if project in [
            "reanalysis-era5-single-levels",
            "reanalysis-era5-single-levels-preliminary-back-extension",
            "reanalysis-era5-pressure-levels",
            "reanalysis-era5-pressure-levels-preliminary-back-extension",
        ]:
            request_kwargs.update(dict(product_type=product))

        if pressure_levels:
            if separate_pressure_level_requests:
                for level in pressure_levels:
                    request_kwargs.update(dict(pressure_level=[level]))
                    netcdf_name = (
                        f"{var}{level}_{timestep}_ecmwf_{'-'.join(project.split('-')[1:])}"
                        f"_{product}_{domain.upper()}_{year}{month}.nc"
                    )

                    if Path(netcdf_name).exists():
                        logging.info(
                            f"Dataset {netcdf_name} already exists. Continuing..."
                        )
                        continue

                    if not dry_run:
                        c = Client()
                        c.retrieve(
                            project,
                            request_kwargs,
                            netcdf_name,
                        )
                    else:
                        logging.info(project)
                        logging.info(request_kwargs)
                        logging.info(netcdf_name)

                continue
            else:
                request_kwargs.update(dict(pressure_level=pressure_levels))

        netcdf_name = (
            f"{var}_{timestep}_ecmwf_{'-'.join(project.split('-')[1:])}"
            f"_{product}_{domain.upper()}_{year}{month}.nc"
        )

        if Path(netcdf_name).exists():
            logging.info(f"Dataset {netcdf_name} already exists. Continuing...")
            continue

        if not dry_run:
            c = Client()
            c.retrieve(
                project,
                request_kwargs,
                netcdf_name,
            )
        else:
            logging.info(project)
            logging.info(request_kwargs)
            logging.info(netcdf_name)


def rename_era5_files(path: Union[os.PathLike, str]) -> None:
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
