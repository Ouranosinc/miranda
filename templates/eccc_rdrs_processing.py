import logging
from pathlib import Path

from miranda.convert.eccc_rdrs import convert_rdrs, rdrs_to_daily
from miranda.io import concat_rechunk_zarr


def main():
    home = Path("~").expanduser()
    dask_dir = home.joinpath("tmpout", "dask")
    dask_dir.mkdir(parents=True, exist_ok=True)
    dask_kwargs = dict(
        n_workers=8,
        threads_per_worker=4,
        memory_limit="7GB",
        dashboard_address=8998,
        local_directory=dask_dir,
        silence_logs=logging.ERROR,
    )

    convert_rdrs(
        project="rdrs-v21",
        input_folder=Path(home).joinpath("RDRS_v21", "caspar"),
        output_folder=Path(home).joinpath("RDRS_v21", "tmp/ECCC/RDRS_v21/NAM"),
        output_format="zarr",
        working_folder=Path(home).joinpath("tmpout", "rdrs"),
        **dask_kwargs,
    )

    rdrs_to_daily(
        project="rdrs-v21",
        input_folder=Path(home).joinpath(
            "RDRS_v21", "tmp/ECCC/RDRS_v21/NAM/1hr"
        ),
        output_folder=Path(home).joinpath("RDRS_v21", "tmp/ECCC/RDRS_v21/NAM/day"),
        working_folder=Path(home).joinpath("tmpout", "rdrs1"),
        overwrite=False,
        year_start=None,
        year_end=None,
        process_variables=None,
        **dask_kwargs,
    )

    # def concat_zarr(infolder=None, outfolder=None, overwrite=False):
    #     list_zarr = sorted(list(infolder.glob("*.zarr")))
    #     outzarr = "_".join(list_zarr[0].stem.split("_")[0:-1])
    #     st_yr = list_zarr[0].stem.split("_")[-1].split("-")[0][0:4]
    #     end_yr = list_zarr[-1].stem.split("_")[-1].split("-")[0][0:4]
    #     outzarr = f"{outzarr}_{st_yr}_{end_yr}.zarr"
    #     outzarr = outfolder.joinpath(outzarr)
    #     print(outzarr)
    #
    #     if not outzarr.exists() or overwrite:
    #         if "day" in infolder.as_posix():
    #             chunks = dict(time=(365 * 4) + 1, rlon=50, rlat=50)
    #         else:
    #             chunks = dict(time=(24 * 30 * 2), rlon=50, rlat=50)
    #
    #         # maketemp files 1 zarr per 4 years
    #         years = [y for y in range(int(st_yr), int(end_yr) + 1)]
    #         years = [years[x: x + 4] for x in range(0, len(years), 4)]
    #         for year in years:
    #             print(year)
    #             list_zarr1 = sorted(
    #                 [
    #                     zarrfile
    #                     for zarrfile in list_zarr
    #                     if int(zarrfile.stem.split("_")[-1].split("-")[0][0:4]) in year
    #                 ]
    #             )
    #             assert len(list_zarr1) / len(year) == 12
    #             ds = xr.open_mfdataset(list_zarr1, parallel=True, engine="zarr")
    #
    #             with Client(**dask_kwargs):
    #                 # if outzarr.exists():
    #                 #     zarr_kwargs = {"append_dim": "time", "consolidated": True}
    #                 # else:
    #                 #     zarr_kwargs = {"consolidated": True}
    #                 tmpzarr = outzarr.parent.joinpath(
    #                     "tmp",
    #                     f"{outzarr.stem.split(f'_{st_yr}_')[0]}_{year[0]}-{year[-1]}.zarr",
    #                 )
    #                 tmpzarr.parent.mkdir(exist_ok=True, parents=True)
    #                 print(f"{year} writing to {tmpzarr.as_posix()}")
    #
    #                 job = delayed_write(
    #                     ds=ds,
    #                     outfile=tmpzarr,
    #                     output_format="zarr",
    #                     target_chunks=chunks,
    #                     overwrite=overwrite,
    #                 )  # kwargs=zarr_kwargs)
    #                 compute(job)
    #
    #         # get tmp zarrs
    #         list_zarr = sorted(list(tmpzarr.parent.glob("*zarr")))
    #         ds = xr.open_mfdataset(list_zarr, engine="zarr")
    #         with Client(**dask_kwargs):
    #             job = delayed_write(
    #                 ds=ds,
    #                 outfile=outzarr,
    #                 output_format="zarr",
    #                 target_chunks=chunks,
    #                 overwrite=overwrite,
    #             )  # kwargs=zarr_kwargs)
    #             compute(job)
    #
    #         shutil.rmtree(tmpzarr.parent)

    # for freq in ["1hr", "day"]:
    #     infolder = Path(home).joinpath("RDRS_v2.1", f"tmp/ECCC/RDRS_v2.1/NAM/{freq}")
    #     for variable in [v for v in infolder.glob("*") if v.is_dir()]:
    #         concat_rechunk_zarr(
    #             input_folder=variable,
    #             output_folder=Path(home).joinpath(
    #                 "RDRS_v2.1", f"converted/ECCC/RDRS_v2.1/NAM/{freq}/{variable.name}"
    #             ),
    #             overwrite=False,
    #         )


if __name__ == "__main__":
    main()
