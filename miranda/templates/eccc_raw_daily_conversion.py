import logging
from functools import partial
from multiprocessing import Pool
from pathlib import Path

from miranda.eccc import aggregate_nc_files
from miranda.eccc import convert_daily_flat_files

if __name__ == "__main__":

    time_step = "daily"
    var_codes = [
        1,
        2,
        3,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    ]
    station_file = "/home/tjs/Desktop/ec_data/Station Inventory EN.csv"
    source_data = Path("/home/tjs/Desktop/ec_data/eccc_all")

    # p = Pool()
    # func = partial(convert_daily_flat_files, source_data, source_data)
    # logging.info(func)
    # p.map(func, var_codes)
    # p.close()
    # p.join()

    # convert_daily_flat_files(
    #     source_files=source_data, output_folder=source_data, variables=var_codes
    # )

    # q = Pool()
    # func = partial(
    #     aggregate_nc_files, source_data, source_data, station_file, time_step
    # )
    # logging.info(func)
    # q.map(func, var_codes)
    # q.close()
    # q.join()

    for var in var_codes:
        # var_name = eccc_cf_daily_metadata(var)["nc_name"]
        # out_file = source_data.joinpath(
        #     "eccc_daily_{}".format(date.today().strftime("%Y%m%d"))
        # )
        aggregate_nc_files(
            source_files=source_data,
            output_folder=source_data,
            variables=var,
            station_inventory=station_file,
            time_step="daily",
        )
