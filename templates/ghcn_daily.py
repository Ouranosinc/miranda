import logging
from pathlib import Path
import os
import shutil
from dask.diagnostics import ProgressBar
from miranda.convert.ghcn import convert_ghcn_bychunks
from miranda.convert.ghcn import chunk_list
import xarray as xr
from tempfile import TemporaryDirectory


def main():

    logging.basicConfig(level=logging.INFO)
    working_folder = Path.home().joinpath("ghcnd")
    start_year = 1981
    end_year = 2020

    # download and convert ghcn data by chunks of nstations
    convert_ghcn_bychunks(
        project="ghcnd",
        working_folder=working_folder,
        lon_bnds=[-76, -74],
        lat_bnds=[44, 46],
        start_year=start_year,
        end_year=end_year,
        overwrite=False,
        delete_raw=False,
        nstations=100,
        n_workers=6,
    )

    # combine zarrs
    for var1 in [v for v in working_folder.joinpath("zarr").iterdir() if v.is_dir()]:
        outfolder = working_folder.joinpath(f"final/{var1.name}")
        outfolder.mkdir(parents=True, exist_ok=True)
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
        for vv in ds.data_vars:
            if "flag" in vv:
                ds[vv] = ds[vv].fillna("").astype("str")

        with ProgressBar():
            ds.chunk(dict(station=250, time=365 * 4 + 1)).to_zarr(
                outzarr.with_suffix(".tmp.zarr"), mode="w"
            )
        shutil.move(outzarr.with_suffix(".tmp.zarr"), outzarr)


if __name__ == "__main__":
    main()
