# import logging

from logging import config

from miranda.scripting import LOGGING_CONFIG

# import numpy as np
# import pandas as pd
# import xarray as xr
# from dask.diagnostics import ProgressBar


config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "aggregate_stations",
    "convert_hourly_csv",
    "convert_daily_csv",
    "merge_converted_variables",
]


def convert_hourly_csv():
    pass


def convert_daily_csv():
    pass


def aggregate_stations():
    pass


def merge_converted_variables():
    pass
