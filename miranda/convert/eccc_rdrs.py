import datetime
import logging
import logging.config
import os
import shutil
import warnings
from pathlib import Path

import xarray as xr
from _data_corrections import file_conversion, load_json_data_mappings
from dask.distributed import Client
from dask import compute
from miranda.convert.utils import daily_aggregation, delayed_write
from miranda.scripting import LOGGING_CONFIG
from miranda.units import get_time_frequency

warnings.simplefilter("ignore")
logging.config.dictConfig(LOGGING_CONFIG)
home = os.path.expanduser("~")
dask_dir = Path(home).joinpath("tmpout", "dask")
dask_dir.mkdir(parents=True, exist_ok=True)
dask_kwargs = dict(
    n_workers=5,
    threads_per_worker=5,
    memory_limit="7GB",
    dashboard_address=8999,
    local_directory=dask_dir,
    silence_logs=logging.ERROR,
)


def main(
    project=None,
    input_folder=None,
    output_folder=None,
    output_format=None,
    working_folder=None,
):
    var_attrs = load_json_data_mappings(project=project)["variable_entry"]
    freq_dict = dict(h="hr", d="day")

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
                ds1 = xr.open_dataset(nc, chunks="auto")
                if ds_allvars is None:
                    ds_allvars = ds1
                    outfreq, meaning = get_time_frequency(ds1)
                    outfreq = (
                        f"{outfreq[0]}{freq_dict[outfreq[1]]}"
                        if meaning == "hour"
                        else freq_dict[outfreq[1]]
                    )

                else:
                    ds_allvars = xr.concat([ds_allvars, ds1], dim="time")
            ds_allvars.attrs["Remarks"] = ds_allvars.attrs["Remarks"].replace(
                "Variable names", "Original variable names"
            )
            for month in range(1, 13):
                ds_month = ds_allvars.sel(time=f"{year}-{str(month).zfill(2)}")
                for vv in var_attrs.keys():
                    drop_vars = _get_drop_vars(ncfile=ncfiles[0], keep_vars=[vv, 'rotated_pole'])

                    # outfile_root = ncfiles[0].stem.split()
                    ds_out = ds_month.drop_vars(drop_vars)

                    overwrite_flag = False

                    write_folder = working_folder if working_folder else output_folder

                    job = file_conversion(
                        input=ds_out,
                        project=project,
                        output_path=write_folder.joinpath(
                            var_attrs[vv]["_cf_variable_name"]
                        ),
                        output_format=output_format,
                        add_version_hashes=False,
                        chunks=dict(time=len(ds_out.time), rlon=50, rlat=50),
                        compute=False,
                        correct_variable_names=True,
                        overwrite=True,
                    )

                    if working_folder:
                        for subdir in [
                            f for f in working_folder.glob("*") if f.is_dir()
                        ]:
                            output_folder.joinpath(outfreq, subdir.name).mkdir(
                                exist_ok=True
                            )
                            for zarr in subdir.glob("*.zarr"):
                                new_zarr = zarr.name.split("_")
                                ## TODO manage output file name
                                new_zarr = f"{new_zarr[0]}_{outfreq}_eccc_{new_zarr[1]}_NAM_{new_zarr[2]}"
                                if (
                                    overwrite_flag
                                    or not output_folder.joinpath(
                                        outfreq, subdir.name, new_zarr
                                    ).exists()
                                ):
                                    with Client(**dask_kwargs):
                                        job.compute()
                                    shutil.copytree(
                                        zarr,
                                        output_folder.joinpath(
                                            outfreq, subdir.name, new_zarr.name
                                        ),
                                        dirs_exist_ok=True,
                                    )
                                else:
                                    print(
                                        output_folder.joinpath(
                                            outfreq, subdir.name, new_zarr
                                        ),
                                        " exists continuing ...",
                                    )
                            shutil.rmtree(subdir)


def get_daily_output_zarrpath(zarrstem, hourly_var, daily_var, year):
    outzarr = zarrstem.replace(hourly_var, daily_var).split("_")
    outzarr[-1] = f"{year}"
    outzarr[1] = "day"
    return f"{'_'.join(outzarr)}.zarr"


def rdrs_to_daily(
    input_folder=None, output_folder=None, working_folder=None, overwrite=False
):
    input_folder

    for var_folder in [v for v in input_folder.glob("*") if v.is_dir()]:
        zarrs = sorted(list(var_folder.glob("*.zarr")))
        year1 = xr.open_zarr(zarrs[0]).time.dt.year.min().values
        year2 = xr.open_zarr(zarrs[-1]).time.dt.year.max().values
        for year in range(year1, year2 + 1):

            infiles = sorted(list(var_folder.glob(f"*_{year}*.zarr")))

            if len(infiles) > 12 or len(infiles) == 0:
                raise ValueError(
                    f"found {len(infiles)} input files ... expected 12 or less"
                )
            outvars = daily_aggregation(xr.open_zarr(infiles[0]), keys_only=True)
            for vv in outvars:
                ds = None
                outzarr = get_daily_output_zarrpath(
                    infiles[0].stem, var_folder.name, vv, year
                )
                outzarr = output_folder.joinpath(vv, outzarr)
                if not outzarr.exists() or overwrite:
                    if ds is None:
                        ds = xr.open_mfdataset(infiles, engine="zarr")
                        ds = daily_aggregation(ds=ds)
                    outzarr.parent.mkdir(parents=True, exist_ok=True)
                    if outzarr.exists():
                        shutil.rmtree(outzarr)

                    if working_folder:
                        tmp_zarr = working_folder.joinpath(outzarr.name)
                        if tmp_zarr.exists():
                            shutil.rmtree(tmp_zarr)
                        with Client(**dask_kwargs):
                            ds[vv].chunk(dict(time=-1, rlon=50, rlat=50)).to_zarr(
                                tmp_zarr
                            )
                        shutil.copytree(tmp_zarr, outzarr, dirs_exist_ok=True)
                        shutil.rmtree(tmp_zarr)
                    else:
                        with Client(**dask_kwargs):
                            ds[var_folder.name].chunk(
                                dict(time=-1, rlon=50, rlat=50)
                            ).to_zarr(outzarr)
                else:
                    print(outzarr.as_posix(), "exists. Continuing..")

def concat_zarr(infolder=None, outfolder=None, overwrite=False):

    list_zarr = sorted(list(infolder.glob("*.zarr")))
    outzarr = "_".join(list_zarr[0].stem.split('_')[0:-1])
    st_yr = list_zarr[0].stem.split('_')[-1].split('-')[0][0:4]
    end_yr = list_zarr[-1].stem.split('_')[-1].split('-')[0][0:4]
    outzarr = f"{outzarr}_{st_yr}_{end_yr}.zarr"
    outzarr = outfolder.joinpath(outzarr)
    print(outzarr)

    if not outzarr.exists() or overwrite:

        if 'day' in infolder.as_posix():
            chunks = dict(time=(365*4)+1, rlon=50, rlat=50)
            chunk_factor = 1
        else:
            chunks = dict(time=(24*30*2), rlon=50, rlat=50)
            chunk_factor = 1

        #maketemp files 1 zarr per 4 years
        years = [y for y in range(int(st_yr), int(end_yr)+1)]
        years = [years[x:x + 4] for x in range(0, len(years), 4)]
        for year in years:
            print(year)
            list_zarr1 = sorted([l for l in list_zarr if int(l.stem.split('_')[-1].split('-')[0][0:4]) in year])
            assert len(list_zarr1) / len(year) == 12
            ds = xr.open_mfdataset(list_zarr1, parallel=True, engine='zarr')

            with Client(**dask_kwargs) as client:
                # if outzarr.exists():
                #     zarr_kwargs = {"append_dim": "time", "consolidated": True}
                # else:
                #     zarr_kwargs = {"consolidated": True}
                tmpzarr = outzarr.parent.joinpath("tmp",f"{outzarr.stem.split(f'_{st_yr}_')[0]}_{year[0]}-{year[-1]}.zarr")
                tmpzarr.parent.mkdir(exist_ok=True, parents=True)
                print(f'{year} writing to {tmpzarr.as_posix()}')

                job = delayed_write(ds=ds, outfile=tmpzarr, output_format="zarr", target_chunks=chunks, overwrite=overwrite)# kwargs=zarr_kwargs)
                compute(job)

        # get tmp zarrs
        list_zarr = sorted(list(tmpzarr.parent.glob('*zarr')))
        ds = xr.open_mfdataset(list_zarr, engine="zarr")
        with Client(**dask_kwargs) as client:
            job = delayed_write(ds=ds, outfile=outzarr, output_format="zarr", target_chunks=chunks,
                                overwrite=overwrite)  # kwargs=zarr_kwargs)
            compute(job)

        shutil.rmtree(tmpzarr.parent)

def _get_drop_vars(ncfile=None, keep_vars=None):
    drop_vars = list(xr.open_dataset(ncfile).data_vars)
    return list(set(drop_vars) - set(keep_vars))


if __name__ == "__main__":
    project = "rdrs-v2.1"
    kwargs = dict(
        project=project,
        input_folder=Path(home).joinpath("RDRS_v2.1", "caspar"),
        output_folder=Path(home).joinpath("RDRS_v2.1", "tmp/ECCC/RDRS_v2.1/NAM"),
        output_format="zarr",
        working_folder=Path(home).joinpath("tmpout", "rdrs"),
    )

    main(**kwargs)

    kwargs = dict(
        input_folder=Path(home).joinpath(
            "RDRS_v2.1", "converted/ECCC/RDRS_v2.1/NAM/1hr"
        ),
        output_folder=Path(home).joinpath(
            "RDRS_v2.1", "tmp/ECCC/RDRS_v2.1/NAM/day"
        ),
        working_folder=Path(home).joinpath("tmpout", "rdrs"),
        overwrite=False,

    )
    rdrs_to_daily(**kwargs)

    for freq in [ '1hr','day',]:
        infolder = Path(home).joinpath("RDRS_v2.1", f"tmp/ECCC/RDRS_v2.1/NAM/{freq}")
        for vv in [i for i in infolder.glob('*') if i.is_dir()]:
            concat_zarr(infolder=vv,
                        outfolder=Path(home).joinpath("RDRS_v2.1", f"converted/ECCC/RDRS_v2.1/NAM/{freq}/{vv.name}"),
                        overwrite=False)
