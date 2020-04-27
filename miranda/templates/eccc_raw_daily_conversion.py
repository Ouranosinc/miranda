import logging
from functools import partial
import itertools as it
from multiprocessing import Pool
from pathlib import Path

from miranda.eccc import aggregate_nc_files
from miranda.eccc import convert_daily_flat_files
from miranda.eccc._raw import _combine_years
if __name__ == "__main__":

    time_step = "daily"
    var_codes = [
        #1,
        #2,
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
    station_file = "/media/sf_VMshare/Trevor/data/Station Inventory EN.csv"
    source_data = Path("/media/sf_VMshare/Trevor/data/netcdf")

    # p = Pool()
    # func = partial(convert_daily_flat_files, source_data, source_data)
    # logging.info(func)
    # p.map(func, var_codes)
    # p.close()
    # p.join()
    #
    # convert_daily_flat_files(
    #     source_files=source_data, output_folder=source_data, variables=var_codes
    # )
    #
    # q = Pool()
    # func = partial(
    #     aggregate_nc_files, source_data, source_data, station_file, time_step
    # )
    # logging.info(func)
    # q.map(func, var_codes)
    # q.close()
    # q.join()


    #TODO add loop on all variables - Do this here or elsewhere??
    # Reduce file number : Combine all years for each station - do for tas as a test
    inrep = '/home/travis/doris_home/logan/scen3/smith/eccc/tas'
    outrep = "/media/sf_VMshare/Trevor/data/netcdf/tas"
    Path(outrep).mkdir(parents=True,exist_ok=True)
    station_dirs= [x for x in Path(inrep).iterdir() if x.is_dir()]
    # stats = [s for s in stats if len(list(Path(outrep).glob(f'{s.name}_*.nc')))==0]
    #
    # combs = list(it.product(*[station_dirs, [outrep]]))
    # for c in combs:
    #     _combine_years(c)
    # q = Pool(16)
    # q.map(_combine_years, combs)
    # q.close()
    # q.join()


    for var in var_codes:
        aggregate_nc_files(
            source_files=source_data,
            output_folder=source_data,
            variables=var,
            station_inventory=station_file,
            time_step="daily",
            mf_dataset_freq='25YS',

        )
