import logging
from pathlib import Path
import shutil
from dask.diagnostics import ProgressBar
from pint.facets.numpy.numpy_func import out_unit

from miranda.convert.stationdata import convert_statdata_bychunks, q_flag_dict
from miranda.eccc._homogenized import download_canhomt
import xarray as xr


def main():

    logging.basicConfig(level=logging.INFO)
    working_folder = Path.home().joinpath("CanHomTV4/CanHomT_dlyV4")

    lon_bnds = [-80, -75]
    lat_bnds = [44, 45]

    n_stations = 100

    start_year = 1800
    end_year = None

    update_raw = True

    # download station data
    download_canhomt(
        project="canhomt_dly",
        working_folder=working_folder,
        update_raw=update_raw,
    )

    # convert ghcn data by chunks of n_stations
    convert_statdata_bychunks(
        project="canhomt_dly",
        working_folder=working_folder,
        lon_bnds=lon_bnds,
        lat_bnds=lat_bnds,
        start_year=start_year,
        end_year=end_year,
        n_stations=n_stations,
        n_workers=6,
    )

    # combine zarrs
    for var1 in [v for v in working_folder.joinpath("zarr").iterdir() if v.is_dir()]:
        out_folder = working_folder.joinpath(f"final/{var1.name}")
        out_folder.mkdir(parents=True, exist_ok=True)
        if out_folder.exists():
            shutil.rmtree(out_folder)
        inzarrs = sorted(var1.glob("*.zarr"))
        ds = xr.concat(
            [xr.open_zarr(z, decode_timedelta=False) for z in inzarrs], dim="station"
        ).sortby(["station", "time"])
        out_zarr = out_folder.joinpath(
            f'{var1.name}_day_GovCan_CamHomT_v4_{ds.time[0].dt.strftime("%Y%m%d").values}-{ds.time[-1].dt.strftime("%Y%m%d").values}.zarr'
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
                out_zarr.with_suffix(".tmp.zarr"), mode="w"
            )
        shutil.move(out_zarr.with_suffix(".tmp.zarr"), out_zarr)


if __name__ == "__main__":
    main()
