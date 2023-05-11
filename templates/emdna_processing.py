import glob
import os
from datetime import datetime
from pathlib import Path

from dask.diagnostics import ProgressBar

from miranda import convert, io, structure
from miranda.decode import Decoder


def preprocess_dna(ds):
    ds = ds.rename_dims({"x": "lon", "y": "lat"})
    ds = ds.set_index(lat="latitude", lon="longitude", time="date")
    ds["time"] = [datetime.strptime(str(d), "%Y%m%d") for d in ds["time"].values]
    # ds= ds.rename({'tmean':'tas', 'trange':'dtr'})
    return ds


def main():
    path = "/jarre/logan/EMDNA/"

    files_by_member = convert.gather_emdna("emdna", path)

    for member, files in files_by_member.items():
        ds = convert.dataset_conversion(
            files, project="EMDNA", preprocess=preprocess_dna
        )
        ds.attrs["member"] = member
    # TODO: deal with units

    for var in ds.data_vars:
        new_path = f"/jarre/scenario/staging/reconstruction/USask/EMDNA_v10/day/{var}/"
        if not os.path.exists(new_path):
            with ProgressBar():
                # write to disk
                io.write_dataset(
                    ds[[var]],
                    project="EMDNA",
                    output_path=new_path,
                    output_format="zarr",
                    overwrite=True,
                    compute=True,
                )
                # TODO: compression ?


if __name__ == "__main__":
    main()
