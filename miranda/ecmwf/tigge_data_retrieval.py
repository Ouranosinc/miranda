# import calendar
import functools
import logging
import multiprocessing
import os
from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path
from typing import List, Optional

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


def prepare_tigge_request(
    *,
    variables: List[str] = None,
    providers: Optional[List[str]] = None,
    forecast_type: str = "pf",
    times: Optional[List[str]] = None,
    dates: Optional[List[str]] = None,
    start_day: Optional[str] = None,
    end_day: Optional[str] = None,
    output_folder: Optional[os.PathLike] = None,
):
    # Providers of interest
    if providers is None:
        providers = ["ecmf"]

    if output_folder is None:
        target = Path().cwd().joinpath("download")
    else:
        target = Path(output_folder)
    Path(target).mkdir(exist_ok=True)
    os.chdir(target)

    if times is None:
        times = ["00/12"]

    if start_day and end_day:
        start = (dt.today() - td(days=4)).strftime("%Y-%m-%d")
        finish = (dt.today() - td(days=3)).strftime("%Y-%m-%d")
        date_range = f"{start}/to/{finish}"
        times.append(date_range)
    elif dates:
        pass
    else:
        raise ValueError()

    tigge_variables = dict()
    tigge_variables["uas"] = 165
    tigge_variables["vas"] = 166
    tigge_variables["tas"] = 167
    tigge_variables["tds"] = 168
    tigge_variables["pr"] = 222
    tigge_variables["swe"] = 228144
    tigge_variables["tasmax"] = 121
    tigge_variables["tasmin"] = 122

    if variables is None:
        variables = tigge_variables.keys()

    project_members = dict(
        dems=44,
        kwbc=26,
        cwao=20,
        ecmf=50,
        babj=30,
        egrr=23,
        rksl=24,
        rjtd=50,
        edzw=40,
        lfpw=34,
        ammc=32,
    )

    for v in variables:
        var_num = tigge_variables[v]
        for t in times:
            for d in dates:
                if forecast_type == "pf":
                    for p in providers:
                        numbers = project_members[p]
                        p = multiprocessing.Pool(processes=4)
                        func = functools.partial(
                            collect_tigge,
                            variable_name=v,
                            variable_code=var_num,
                            time=t,
                            provider=p,
                            numbers=numbers,
                            dates=d,
                        )

                        logging.info([func, dt.now().strftime("%Y-%m-%d %X")])

                        p.map(func, dates)
                        p.close()
                        p.join()

                elif forecast_type == "cf":
                    p = multiprocessing.Pool(processes=4)
                    func = functools.partial(
                        collect_tigge,
                        variable_name=v,
                        variable_code=var_num,
                        time=t,
                        provider=p,
                        dates=d,
                    )

                    logging.info([func, dt.now().strftime("%Y-%m-%d %X")])

                    p.map(func, dates)
                    p.close()
                    p.join()


def collect_tigge(
    variable_name: str,
    variable_code: str,
    time: str,
    forecast_type: str,
    provider: str,
    numbers: Optional[int],
    date: str,
):
    numbers = "/".join([str(n) for n in range(1, numbers + 1)])
    output_name = (
        f"{variable_name}_{provider}_{'-'.join(time.split('/'))}_tigge_reanalysis_6h_"
        f"{date.split('/')[0]}_{date.split('/')[-1]}.grib2"
    )

    # Note: This is only valid for ECMWF at the moment.
    steps = (
        "0/6/12/18/24/30/36/42/48/54/60/66/72/78/84/90/96/102/108/114/120/126/132/138/144/150/156/162/168/"
        "174/180/186/192/198/204/210/216/222/228/234/240/246/252/258/264/270/276/282/288/294/300/306/312/318/324/"
        "330/336/342/348/354/360"
    )

    # Remove steps 0 for tasmax and tasmin
    if variable_code in [121, 122]:
        steps = steps[2:]

    request = {
        "class": "ti",
        "dataset": "tigge",
        "date": date,
        "expver": "prod",
        "grid": "0.5/0.5",
        "levtype": "sfc",
        "origin": provider,
        "param": variable_code,
        "step": steps,
        "time": time,
        "type": forecast_type,
        "target": output_name,
    }
    if numbers:
        request.update({"number": numbers})

    server = ECMWFDataServer()
    server.retrieve(request)
