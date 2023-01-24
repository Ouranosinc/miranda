from os import getenv
from pathlib import Path

from miranda.structure import structure_datasets

if __name__ == "__main__":
    in_files = getenv("in")
    out_files = getenv("out")

    input_path = Path(in_files)
    output_path = Path(out_files)

    structure_datasets(
        input_path,
        output_path,
        project="converted",
        guess=False,
        method="copy",
        make_dirs=True,
        filename_pattern="*.zarr",
    )
