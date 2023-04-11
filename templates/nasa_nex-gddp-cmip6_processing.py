from pathlib import Path

from miranda import convert, io


def main():
    path_nex_out = Path("/jarre/scenario/jlavoie/data/NEX-GDDP-CMIP6_downloaded/")
    nex_files = convert.gather_nex(path_nex_out)

    for path, list_files in nex_files.items():
        ds = convert.dataset_conversion(
            list_files,
            add_version_hashes=False,
            preprocess=False,
            project="NEX-GDDP-CMIP6",
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
