
import logging.config
import logging
import os
import shutil
from pathlib import Path
import datetime
import xarray as xr
from dask.distributed import Client
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency
from miranda.convert.utils import daily_aggregation
from _data_corrections import load_json_data_mappings, file_conversion
import warnings
warnings.simplefilter('ignore')
logging.config.dictConfig(LOGGING_CONFIG)
home = os.path.expanduser('~')
dask_dir = Path(home).joinpath('tmpout', "dask")
dask_dir.mkdir(parents=True, exist_ok=True)
dask_kwargs = dict(n_workers=10, threads_per_worker=3, memory_limit="3GB", dashboard_address=8999,
                   local_directory=dask_dir, silence_logs=logging.ERROR)

def main(project=None, input_folder=None, output_folder=None, output_format=None, working_folder=None):
    var_attrs = load_json_data_mappings(project=project)['variable_entry']
    freq_dict= dict(h="hr", d="day")

    if working_folder:
        if working_folder.exists():
            shutil.rmtree(working_folder)
        working_folder.mkdir(parents=True)

    year = datetime.date.today().year
    for year in range(1950, year):
        if not list(input_folder.glob(f"{year}*.nc")):
            print(f"no files found for year {year} ... continuing")
            continue
        # Time stamps starts at noon and flow into subsequent months!!
        # Need full year plus previous december in order to easily produce complete hourly frequency monthly files
        ncfiles = sorted(list(input_folder.glob(f"{year - 1}12*.nc")))
        if ncfiles:
            ncfiles = [ncfiles[-1]]
        ncfiles.extend(sorted(list(input_folder.glob(f"{year}*.nc"))))
        ds_allvars = None
        if ncfiles:
            for nc in ncfiles:
                print(nc.name)
                ds1 = xr.open_dataset(nc, chunks='auto')
                if ds_allvars is None:
                    ds_allvars = ds1
                    outfreq, meaning = get_time_frequency(ds1)
                    outfreq = f"{outfreq[0]}{freq_dict[outfreq[1]]}" if meaning == 'hour' else freq_dict[outfreq[1]]

                else:
                    ds_allvars = xr.concat([ds_allvars, ds1], dim='time')
            ds_allvars.attrs['Remarks'] = ds_allvars.attrs['Remarks'].replace('Variable names', 'Original variable names')
            for month in range(1, 13):
                ds_month = ds_allvars.sel(time=f"{year}-{str(month).zfill(2)}")
                for vv in var_attrs.keys():
                    drop_vars = _get_drop_vars(ncfile=ncfiles[0], keep_vars=[vv])

                    # outfile_root = ncfiles[0].stem.split()
                    ds_out = ds_month.drop_vars(drop_vars)

                    overwrite_flag = False

                    write_folder = working_folder if working_folder else output_folder

                    job = file_conversion(input=ds_out, project=project,
                                    output_path=write_folder.joinpath(var_attrs[vv]['_cf_variable_name']),
                                    output_format=output_format,
                                    add_version_hashes=False, chunks=dict(time=len(ds_out.time), rlon=50, rlat=50),
                                    compute=False, correct_variable_names=True, overwrite=True)

                    if working_folder:
                        for subdir in [f for f in working_folder.glob('*') if f.is_dir()]:
                            output_folder.joinpath(outfreq, subdir.name).mkdir(exist_ok=True)
                            for zarr in subdir.glob('*.zarr'):
                                new_zarr = zarr.name.split('_')
                                ## TODO manage output file name
                                new_zarr = f"{new_zarr[0]}_{outfreq}_eccc_{new_zarr[1]}_NAM_{new_zarr[2]}"
                                if overwrite_flag or not output_folder.joinpath(outfreq, subdir.name, new_zarr).exists():
                                    with Client(**dask_kwargs):
                                        job.compute()
                                    shutil.copytree(zarr, output_folder.joinpath(outfreq, subdir.name, new_zarr.name), dirs_exist_ok=True)
                                else:
                                    print(output_folder.joinpath(outfreq, subdir.name, new_zarr), ' exists continuing ...')
                            shutil.rmtree(subdir)
def get_daily_output_zarrpath(zarrstem, hourly_var, daily_var, year):
    outzarr = zarrstem.replace(hourly_var, daily_var).split('_')
    outzarr[-1] = f"{year}"
    outzarr[1] = "day"
    return f"{'_'.join(outzarr)}.zarr"
def rdrs_to_daily(input_folder=None, output_folder=None, working_folder=None, overwrite=False):
    input_folder

    for var_folder in [v for v in input_folder.glob("*") if v.is_dir()]:
        zarrs = sorted(list(var_folder.glob("*.zarr")))
        year1 = xr.open_zarr(zarrs[0]).time.dt.year.min().values
        year2 = xr.open_zarr(zarrs[-1]).time.dt.year.max().values
        for year in range(year1, year2+1):

            infiles = sorted(list(var_folder.glob(f"*_{year}*.zarr")))

            if len(infiles)>12 or len(infiles)==0:
                raise ValueError(f"found {len(infiles)} input files ... expected 12 or less")
            outvars = daily_aggregation(xr.open_zarr(infiles[0]), keys_only=True)
            for vv in outvars:
                ds = None
                outzarr = get_daily_output_zarrpath(infiles[0].stem, var_folder.name, vv, year)
                outzarr = output_folder.joinpath(vv, outzarr)
                if not outzarr.exists() or overwrite:
                    if ds is None:
                        ds = xr.open_mfdataset(infiles, engine='zarr')
                        ds = daily_aggregation(ds=ds)
                    outzarr.parent.mkdir(parents=True, exist_ok=True)
                    if outzarr.exists():
                        shutil.rmtree(outzarr)

                    if working_folder:
                        tmp_zarr = working_folder.joinpath(outzarr.name)
                        if tmp_zarr.exists():
                            shutil.rmtree(tmp_zarr)
                        with Client(**dask_kwargs):
                            ds[vv].chunk(dict(time=-1, rlon=50, rlat=50)).to_zarr(tmp_zarr)
                        shutil.copytree(tmp_zarr, outzarr, dirs_exist_ok=True)
                        shutil.rmtree(tmp_zarr)
                    else:
                        with Client(**dask_kwargs):
                            ds[var_folder.name].chunk(dict(time=-1, rlon=50, rlat=50)).to_zarr(outzarr)
                else:
                    print(outzarr.as_posix(), "exists. Continuing..")


def _get_drop_vars(ncfile=None, keep_vars = None):
    drop_vars = list(xr.open_dataset(ncfile).data_vars)
    return list(set(drop_vars) - set(keep_vars))



if __name__ == '__main__':
    project = "rdrs-v2.1"
    kwargs = dict(project=project,
                  input_folder=Path(home).joinpath('RDRS_v2.1', 'caspar'),
                  output_folder=Path(home).joinpath('RDRS_v2.1', 'converted/ECCC/RDRS_v2.1/NAM'),
                  output_format="zarr",
                  working_folder=Path(home).joinpath('tmpout','rdrs')
                  )

    #main(**kwargs)
    rdrs_to_daily(
        input_folder=Path(home).joinpath('RDRS_v2.1', 'converted/ECCC/RDRS_v2.1/NAM/1hr'),
        output_folder=Path(home).joinpath('RDRS_v2.1', 'converted/ECCC/RDRS_v2.1/NAM/day'),
        working_folder=Path(home).joinpath('tmpout', 'rdrs'),
        overwrite=False
    )
