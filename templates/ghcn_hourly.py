import logging
from pathlib import Path
import os
import shutil
from dask.diagnostics import ProgressBar
from miranda.convert.ghcn import convert_ghcn_bychunks, download_ghcn, q_flag_dict
from miranda.convert.ghcn import chunk_list
import xarray as xr
from tempfile import TemporaryDirectory


def main():

    logging.basicConfig(level=logging.INFO)
    working_folder = Path.home().joinpath("ghcnh")
    start_year = 2001
    end_year = 2020

    lon_bnds = [-76, -74]
    lat_bnds = [44, 46]

    nstations = 100
    update_raw = True

    # download station data
    # download_ghcn(
    #     project="ghcnh",
    #     working_folder=working_folder,
    #     lon_bnds=lon_bnds,
    #     lat_bnds=lat_bnds,
    #     update_raw=update_raw,
    # )

    # convert ghcn data by chunks of nstations
    convert_ghcn_bychunks(
        project="ghcnh",
        working_folder=working_folder,
        start_year=start_year,
        end_year=end_year,
        update_from_raw=update_raw,
        nstations=nstations,
        n_workers=6,
    )

    # combine zarrs
    for var1 in [v for v in working_folder.joinpath("zarr").iterdir() if v.is_dir()]:
        outfolder = working_folder.joinpath(f"final/{var1.name}")
        outfolder.mkdir(parents=True, exist_ok=True)
        if outfolder.exists():
            shutil.rmtree(outfolder)
        inzarrs = sorted(var1.glob("*.zarr"))
        ds = xr.concat(
            [xr.open_zarr(z, decode_timedelta=False) for z in inzarrs], dim="station"
        ).sortby(["station", "time"])
        outzarr = outfolder.joinpath(
            f'{var1.name}_day_NOAA_ghcnd_{ds.time[0].dt.strftime("%Y%m%d").values}-{ds.time[-1].dt.strftime("%Y%m%d").values}.zarr'
        )
        # remove encoding
        for c in ds.coords:
            ds[c].load()
            ds[c].encoding = {}
        for c in ds.data_vars:
            ds[c].encoding = {}
            if "flag" in c:
                mask = ds[c].isin(list(q_flag_dict.keys()))
                ds[c] = ds[c].where(mask, "")
                ds[c] = ds[c].astype(str)
        with ProgressBar():
            ds.chunk(dict(station=250, time=365 * 4 + 1)).to_zarr(
                outzarr.with_suffix(".tmp.zarr"), mode="w"
            )
        shutil.move(outzarr.with_suffix(".tmp.zarr"), outzarr)


if __name__ == "__main__":
    main()
