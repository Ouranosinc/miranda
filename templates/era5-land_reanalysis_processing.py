from pathlib import Path

from miranda import convert

path_era5_land_out = Path("~/Desktop").expanduser()
era5_land_files = convert.gather_ecmwf("era5-land", path_era5_land_out)

ERA5_VARIABLES = [
    "d2m",
    "pev",
    "rsn",
    "sde",
    "sd",
    "sf",
    "t2m",
    "tp",
    "u10",
    "v10",
]

convert.reanalysis_processing(
    era5_land_files,
    path_era5_land_out,
    variables=ERA5_VARIABLES,
    aggregate=False,
    target_chunks=None,
    domains="NAM",
    output_format="zarr",
    overwrite=False,
)
