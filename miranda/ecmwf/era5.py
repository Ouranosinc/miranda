import functools
import logging.config
import multiprocessing
import os
from datetime import date
from datetime import datetime as dt
from pathlib import Path
from typing import List, Mapping, Optional, Tuple

from cdsapi import Client

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["request_era5"]


def request_era5(
    variables: Optional[Mapping[str, str]] = None,
    projects: List[str] = None,
    domain: str = "AMNO",
    year_start: int = 1981,
    year_end: Optional[int] = None,
    processes: int = 4,
) -> None:
    """Request ERA5/ERA5-Land from Copernicus Data Store in NetCDF4 format.

    Parameters
    ----------
    variables: Mapping[str, str], optional
    projects : List[{"era5", "era5-land"}]
    domain : {"GLOBAL", "AMNO", "CAN", "QC"}
    year_start : int
    year_end : int, optional
    processes : int

    Returns
    -------
    None
    """
    # Variables of interest
    variable_reference = dict(
        pr="total_precipitation",
        vas="10m_v_component_of_wind",
        uas="10m_u_component_of_wind",
        td="2m_dewpoint_temperature",
        tas="2m_temperature",
        potevap="potential evaporation",
        snd="snow_depth",
        prsn="snowfall",
    )

    v_requested = dict()
    if variables:
        for v in variables:
            if v in variable_reference:
                v_requested[v] = variable_reference[v]
    else:
        v_requested = variable_reference

    if year_end is None:
        year_end = date.today().year
    years = range(year_start, year_end)

    months = [str(d).zfill(2) for d in range(13)]
    yearmonth = list()
    for y in years:
        for m in months:
            yearmonth.append((y, m))

    project_names = list()
    if "era5" in projects:
        project_names.append("reanalysis-era5-single-levels")
    if "era5-land" in projects:
        project_names.append("reanalysis-era5-land")
    product = project_names[0].split("-")[0]

    target = Path().cwd().joinpath("downloaded")

    Path(target).mkdir(exist_ok=True)
    os.chdir(target)

    for p in projects:
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
    year, month = yearmonth
    days = [str(d).zfill(2) for d in range(32)]
    times = ["{}:00".format(str(t).zfill(2)) for t in range(24)]

    if domain.upper() == "GLOBAL":
        region = [90, -180, -90, 180]
    elif domain.upper() == "AMNO":
        domain = "NAM"
        region = [90, -180, 10, -10]
    elif domain.upper() == "CAN":
        region = [83.5, -141, 41.5, -52.5]
    elif domain.upper() == "QC":
        region = [63, -80, 44.5, -57]
    else:
        raise ValueError()

    c = Client()

    for var in variables.keys():
        netcdf_name = f"{var}_{'-'.join(project.split('-')[1:])}_{product}_hourly_{year}{month}_{domain.upper()}.nc"

        if Path(netcdf_name).exists():
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
