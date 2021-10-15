# import calendar
import functools
import logging
import multiprocessing
import os
from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path

from ecmwfapi import ECMWFDataServer

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


def main(variables, members: dict, start_times, dates):

    target = Path().cwd().joinpath("download")

    Path(target).mkdir(exist_ok=True)
    os.chdir(target)

    for variable_name, variable_code in variables.items():
        for member, number in members.items():
            for time in start_times:
                p = multiprocessing.Pool(processes=4)
                func = functools.partial(
                    collect_tigge, variable_name, variable_code, time, member, number
                )

                logging.info([func, dt.now().strftime("%Y-%m-%d %X")])

                p.map(func, dates)
                p.close()
                p.join()


def collect_tigge(
    variable_name: str,
    variable_code: str,
    time: str,
    member: str,
    number: int,
    date: str,
):
    numbers = "/".join([str(n) for n in range(1, number + 1)])
    output_name = (
        f"{variable_name}_{member}_{'-'.join(time.split('/'))}_tigge_reanalysis_6h_"
        f"{date.split('/')[0]}_{date.split('/')[-1]}.grib2"
    )

    steps = (
        "0/6/12/18/24/30/36/42/48/54/60/66/72/78/84/90/96/102/108/114/120/126/132/138/144/150/156/162/168/"
        "174/180/186/192/198/204/210/216/222/228/234/240/246/252/258/264/270/276/282/288/294/300/306/312/318/324/"
        "330/336/342/348/354/360"
    )

    # Remove steps 0 for tasmax and tasmin
    if variable_code in [121, 122]:
        steps = steps[2:]

    server = ECMWFDataServer()
    server.retrieve(
        {
            "class": "ti",
            "dataset": "tigge",
            "date": date,
            "expver": "prod",
            "grid": "0.5/0.5",
            "levtype": "sfc",
            "number": numbers,
            "origin": member,
            "param": variable_code,
            "step": steps,
            "time": time,
            "type": "pf",
            "target": output_name,
        }
    )


if __name__ == "__main__":
    # Variables of interest
    variables = dict()
    variables["uas"] = 165
    variables["vas"] = 166
    variables["tas"] = 167
    variables["tds"] = 168
    variables["pr"] = 222
    variables["swe"] = 228144
    variables["tasmax"] = 121
    variables["tasmin"] = 122

    # Sources and number of members
    members = dict(
        # dems=44,
        # kwbc=26,
        # cwao=20,
        ecmf=50,
        # babj=30,
        # egrr=23,
        # rksl=24,
        # rjtd=50,
        # edzw=40,
        # lfpw=34,
        # ammc=32,
    )

    # Forecast start dates
    start_times = ["00/12"]

    dates = list()
    # Grab entire months
    # years = range(2020, 2021)
    # for year in years:
    #     for month in range(1, 13):
    #         dates.append(
    #             f"{year}-{str(month).zfill(2)}-01/to/{year}-{str(month).zfill(2)}-"
    #             f"{calendar.monthrange(year, month)[1]}"
    #         )

    # Collect individual days
    start = (dt.today() - td(days=4)).strftime("%Y-%m-%d")
    finish = (dt.today() - td(days=3)).strftime("%Y-%m-%d")
    date_range = f"{start}/to/{finish}"
    dates.append(date_range)

    main(variables, members, start_times, dates)
