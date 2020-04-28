from pathlib import Path
from miranda.eccc import aggregate_nc_files
from pathlib import Path

from miranda.eccc import aggregate_nc_files

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
    # station_file = "/home/tjs/Desktop/ec_data/Station Inventory EN.csv"
    station_file = "/media/sf_VMshare/Trevor/data/Station Inventory EN.csv"
    source_data = Path("/home/travis/doris_home/logan/scen3/smith/eccc")

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
    Path(outrep).mkdir(parents=True, exist_ok=True)

    station_dirs = [x for x in Path(source_data).joinpath(variable).iterdir() if x.is_dir()]

    # stats = [s for s in stats if len(list(Path(outrep).glob(f'{s.name}_*.nc')))==0]
    #
    # combs = list(it.product(*[[variable],station_dirs, [outrep]]))
    #
    # for c in combs:
    #     _combine_years(c)
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
