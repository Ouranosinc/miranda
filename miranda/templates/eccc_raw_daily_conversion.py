import logging
from datetime import date
from functools import partial
from multiprocessing import Pool
from pathlib import Path

from miranda.eccc import aggregate_nc_files
from miranda.eccc import convert_daily_flat_files
from miranda.utils import eccc_cf_daily_metadata

if __name__ == "__main__":

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

    func = partial(convert_daily_flat_files, source_data, source_data)
    logging.info(func)

    p = Pool()
    p.map(func, var_codes)
    p.close()
    p.join()

    # convert_daily_flat_files(
    #     source_files=source_data, output_folder=source_data, variables=var_codes
    # )

    for var in var_codes:
        var_name = eccc_cf_daily_metadata(var)["nc_name"]
        out_file = source_data.joinpath(
            "eccc_daily_{}".format(date.today().strftime("%Y%m%d"))
        )
        aggregate_nc_files(
            source_files=source_data,
            output_file=out_file,
            variables=var,
            station_inventory=station_file,
            time_step="daily",
        )
