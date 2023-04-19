from pathlib import Path

from miranda import convert, io


def main():
    path_era5_land_out = Path("~/Desktop").expanduser()
    era5_land_files = convert.gather_ecmwf("era5-land", path_era5_land_out)

    ds = convert.dataset_conversion(
        era5_land_files,
        project="era5-land-monthly-means",
    )

    io.write_dataset(
        ds,
        project="era5-land-monthly-means",
        output_path="~/Desktop/",
        output_format="netcdf",
        overwrite=True,
        compute=True,
    )


if __name__ == "__main__":
    main()
