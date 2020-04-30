from pathlib import Path
import itertools

from miranda.eccc import aggregate_nc_files
from miranda.eccc._raw import _combine_years

if __name__ == "__main__":

    time_step = "hourly"
    var_codes = [
        76,
        77,
        78,
        79,
        80,
        89,
        94,
        123,
        262,
        263,
        264,
        265,
        266,
        267,
        268,
        269,
        270,
        271,
        272,
        273,
        274,
        275,
        276,
        277,
        278,
        279,
        280,
    ]
    # station_file = "/media/sf_VMshare/Trevor/data/Station Inventory EN.csv"
    # source_data = Path("/home/travis/doris_home/logan/scen3/smith/eccc")
    source_data = Path("/home/tjs/Desktop/ec/")
    station_file = source_data.joinpath("Station Inventory EN.csv")

    # p = Pool()
    # func = partial(convert_hourly_flat_files, source_data, source_data)
    # logging.info(func)
    # p.map(func, var_codes)
    # p.close()
    # p.join()
    #
    # convert_hourly_flat_files(
    #    source_files=source_data, output_folder=source_data, variables=var_codes
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

    variable = 'precipitation'
    outrep = Path("/media/sf_VMshare/Trevor/data/netcdf").joinpath(variable)
    outrep.mkdir(parents=True, exist_ok=True)

    station_dirs = [x for x in source_data.joinpath(variable).iterdir() if x.is_dir()]

    combs = list(itertools.product(*[[variable], station_dirs, [outrep]]))
    for c in combs:
        _combine_years(c)
    # q = Pool(16)
    # q.map(_combine_years,combs)
    # q.close()
    # q.join()

    source_data = outrep.parent

    for var in var_codes:
        aggregate_nc_files(
            source_files=source_data,
            output_folder=source_data,
            variables=var,
            station_inventory=station_file,
        )