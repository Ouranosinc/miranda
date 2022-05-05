import logging
from datetime import datetime as dt
from os import getenv
from pathlib import Path

from miranda.ecmwf import request_era5

# logging file configuration
logging.basicConfig(
    filename="{}_{}.log".format(dt.now().strftime("%Y%m%d"), Path(__file__).stem),
    level=logging.INFO,
    datefmt="%H:%M:%S",
)


def main():
    out_files = getenv("out")

    target_folder = Path(out_files).expanduser()

    variables = {
        "d2m": "2m_dewpoint_temperature",
        "pev": "potential_evaporation",
        "rsn": "snow_density",
        "sd": "snow_depth",
        "sf": "snowfall",
        "t2m": "2m_temperature",
        "tp": "total_precipitation",
        "u10": "10m_u_component_of_wind",
        "v10": "10m_v_component_of_wind",
    }
    projects = [
        "era5-single-levels",
        "era5-single-levels-preliminary-back-extension",
    ]

    request_era5(projects, variables=variables, output_folder=target_folder)


if __name__ == "__main__":
    main()