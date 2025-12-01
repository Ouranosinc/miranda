from pathlib import Path
import shutil
import xarray as xr
from dask.diagnostics import ProgressBar

from miranda.convert.NRCanMET import convert_NRCanMET
from miranda.convert._data_corrections import load_json_data_mappings

def main():
    infolder = Path.home().joinpath("data/NRCanMET/netcdf")
    outfolder = Path.home().joinpath("data/miranda_output/NRCanMET/")
    project = "NRCanMET"
    for var in load_json_data_mappings(project=project)['variables']:
        print(f"Variable to convert: {var}")
        # only convert 1 year for template example
        for nc in sorted(list(infolder.glob(f"{var}*_2001_*.nc"))):
            print(f"Processing file: {nc.name}")
            ds = convert_NRCanMET(
                infile= nc,
            )

            vv = list(ds.data_vars)[0]
            out_zarr = outfolder.joinpath(
                f'{vv}_day_GovCan_{project}_{ds.time[0].dt.strftime("%Y%m%d").values}-{ds.time[-1].dt.strftime("%Y%m%d").values}.zarr'
            )

            # Ensure output directory exists
            out_zarr.parent.mkdir(parents=True, exist_ok=True)

            # remove encoding to avoid issues
            for c in ds.coords:
                ds[c].load()
                ds[c].encoding = {}
            for c in ds.data_vars:
                ds[c].encoding = {}

            # Remove calendar attribute from time coordinate to prevent encoding conflicts
            if "calendar" in ds.time.attrs:
                del ds.time.attrs["calendar"]

            with ProgressBar():
                ds.chunk(dict(time=(365*4) + 1, lon=50, lat=50)).to_zarr(
                    out_zarr.with_suffix(".tmp.zarr"), mode="w", zarr_format=2
                )
            shutil.move(out_zarr.with_suffix(".tmp.zarr"), out_zarr)

if __name__ == "__main__":
    main()

