import functools
import logging
import multiprocessing
import os
from datetime import datetime as dt
from pathlib import Path

from cdsapi import Client

# logging file configuration
logging.basicConfig(
    filename="{}_{}.log".format(dt.now().strftime("%Y%m%d"), Path(__file__).stem),
    level=logging.INFO,
    datefmt="%H:%M:%S",
)

# set up logging to console
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# format output to the console
formatter = logging.Formatter("%(name)s : %(asctime)s :  %(levelname)s : %(message)s")
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger("").addHandler(console)
logger = logging.getLogger(__name__)


def main(variables, years):

    months = [str(d).zfill(2) for d in range(13)]
    yearmonth = list()
    for y in years:
        for m in months:
            yearmonth.append((y, m))

    project = "reanalysis-era5-land"
    # project = "reanalysis-era5-single-levels"
    product = project.split("-")[0]

    target = Path().cwd().joinpath("downloaded")

    Path(target).mkdir(exist_ok=True)
    os.chdir(target)

    p = multiprocessing.Pool(processes=4)
    func = functools.partial(collect_era, variables, project, product)

    logging.info([func, dt.now().strftime("%Y-%m-%d %X")])

    p.map(func, yearmonth)
    p.close()
    p.join()


def collect_era(variables, project, product, yearmonth):
    year, month = yearmonth
    days = [str(d).zfill(2) for d in range(32)]
    times = ["{}:00".format(str(t).zfill(2)) for t in range(24)]

    region = "90/-180/10/-10"

    c = Client()

    for var in variables.keys():
        netcdf_name = f"{var}_{'-'.join(project.split('-')[1:])}_{product}_hourly_{year}{month}_NAM.nc"
        if os.path.exists(netcdf_name):
            continue

        c.retrieve(
            project,
            {
                # "product_type": product,
                "variable": variables[var],
                "year": year,
                "month": month,
                "day": days,
                "time": times,
                "area": region,
                "format": "netcdf",
            },
            netcdf_name,
        )


if __name__ == "__main__":
    # Variables of interest
    variables = {
        "pr": "total_precipitation",
        "vas": "10m_v_component_of_wind",
        "uas": "10m_u_component_of_wind",
        "td": "2m_dewpoint_temperature",
        "tas": "2m_temperature",
        "potevap": "potential evaporation",
        "snd": "snow_depth",
        "prsn": "snowfall",
    }

    years = range(1981, 2021)

    main(variables, years)
