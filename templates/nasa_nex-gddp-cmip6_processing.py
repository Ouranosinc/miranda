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
            preprocess=False,
            project="NEX-GDDP-CMIP6",
        )
        # the path is already ok, so just reuse it
        end_path = path.split("downloaded/")[1]

        # get right chunks
        calendar, size = io._rechunk.get_time_attrs(ds)
        target_chunks = io._rechunk.translate_time_chunk(
            {"time": "4 years", "lat": 50, "lon": 50}, calendar, size
        )

        # write to disk
        io.write_dataset(
            ds,
            project="NEX-GDDP-CMIP6",
            output_path=f"OUTPUT_PATH/{end_path}",
            output_format="zarr",
            overwrite=True,
            compute=True,
            chunks=target_chunks,
        )


if __name__ == "__main__":
    main()
