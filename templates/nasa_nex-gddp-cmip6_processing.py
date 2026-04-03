from pathlib import Path

from miranda import convert, io


def main():
    path_nex_out = Path("INPUT_PATH")
    # put all nc that should be in 1 zarr together
    nex_files = convert.gather_nex(path_nex_out)

    for path, list_files in nex_files.items():
        # open as dataset
        ds = convert.dataset_conversion(
            list_files,
            add_version_hashes=False,
            project="NEX-GDDP-CMIP6",
        )
        # the path is already ok, so just reuse it
        end_path = path.split("downloaded/")[1]

        # write to disk
        io.write_dataset(
            ds,
            project="NEX-GDDP-CMIP6",
            output_path=f"OUTPUT_PATH/{end_path}",
            output_format="zarr",
            overwrite=True,
            compute=True,
        )


if __name__ == "__main__":
    main()
