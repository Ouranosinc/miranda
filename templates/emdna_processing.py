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
    return ds


def main():
    path = "/jarre/logan/EMDNA/"

    files_by_member = convert.gather_emdna(path)
    for member, files in files_by_member.items():
        ds = convert.dataset_conversion(
            files, project="EMDNA", preprocess=preprocess_dna
        )
        ds.attrs["member"] = member
        # TODO: deal with units
        facets = ds.attrs
        for var in ds.data_vars:
            facets["variable"] = var
            out_path = "/jarre/scenario/staging/"
            out_path = "/jarre/scenario/jlavoie/test"

            # get full path
            new_path = structure.build_path_from_schema(facets, out_path)
            # add version by hand for now
            new_path = Path(str(new_path).replace("EMDNA", "EMDNA_v10"))

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


if __name__ == "__main__":
    main()
