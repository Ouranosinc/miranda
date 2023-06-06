from __future__ import annotations

import functools
import logging.config
import multiprocessing
import os
from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["request_tigge"]


def request_tigge(
    variables: list[str],
    providers: list[str] | None = None,
    *,
    forecast_type: str = "pf",
    times: list[str] | None = None,
    dates: list[str] | None = None,
    date_start: str | None = None,
    date_end: str | None = None,
    output_folder: os.PathLike | None = None,
    processes: int = 4,
) -> None:
    """Request tigge data from ECMWF in grib format.

    Parameters
    ----------
    variables : list of str
    providers : ist of str, optional
    forecast_type: {"pf", "cf"}
    times : list of str, optional
    dates : list of str. optional
    date_start : str, optional
    date_end : str, optional
    output_folder : os.PathLike, optional
    processes : int

    Returns
    -------
    None
    """

    def _request_direct_tigge(
        variable_name: str,
        variable_code: str,
        time: str,
        fc_type: str,
        provider: str,
        nums: int | None,
        date: str,
    ):
        """Launch formatted request."""
        try:
            from ecmwfapi import ECMWFDataServer
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                f"{_request_direct_tigge.__name__} requires additional dependencies. "
                "Please install them with `pip install miranda[remote]`."
            )

        number_range = ""
        if nums:
            number_range = "/".join([str(n) for n in range(1, nums + 1)])
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
            "type": fc_type,
            "target": output_name,
        }
        if nums:
            request.update({"number": number_range})

        server = ECMWFDataServer()
        server.retrieve(request)

    # Providers of interest
    if providers is None:
        providers = ["ecmf"]

    if output_folder is None:
        target = Path().cwd().joinpath("download")
    else:
        target = Path(output_folder)
    Path(target).mkdir(exist_ok=True)
    os.chdir(target)

    if not times:
        times = ["00/12"]

    if date_start and date_end:
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
                for p in providers:
                    proc = multiprocessing.Pool(processes=processes)
                    config = dict(
                        variable_name=v,
                        variable_code=var_num,
                        time=t,
                        provider=p,
                        dates=d,
                    )
                    num_dict = dict()
                    if forecast_type == "pf":
                        numbers = project_members[p]
                        num_dict[numbers] = numbers

                    func = functools.partial(
                        _request_direct_tigge, **config, **num_dict
                    )

                    logging.info([func, dt.now().strftime("%Y-%m-%d %X")])

                    proc.map(func, dates)
                    proc.close()
                    proc.join()
