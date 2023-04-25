import glob
import os
from pathlib import Path

from miranda import convert, io, structure
from miranda.decode import Decoder


def main():
    # paths = {'ESPO-G6-R2': ["/jarre/scenario/jlavoie/ESPO-G6/FINAL/NAM-rdrs/*",
    #         "/jarre/braun/data/ESPO-G6/final/*"],
    #         'ESPO-G6-E5L': ['/crue/jlavoie/info-crue-cmip6/FINAL_QC/*', 'path_sarah']}

    paths = {"ESPO-G6-R2": ["PATH1", "PATH2"], "ESPO-G6-E5L": ["PATH1", "PATH2"]}

    for project in paths:
        d = Decoder(project)

        files = []
        for p in paths[project]:
            files.extend(glob.glob(p))

        for f in files:
            if project == "ESPO-G6-R2":
                facets = d.decode_espo_g6_r2(Path(f))  # can't make just decode work
            elif project == "ESPO-G6-E5l":
                facets = d.decode_espo_g6_e5l(Path(f))  # can't make just decode work
            else:
                print("Wrong project!")
                break
            new_path = structure.build_path_from_schema(facets, "OUTPUT_PATH/")

            if not os.path.exists(new_path):  # and path not in skip:
                # open as dataset
                ds = convert.dataset_conversion(
                    f,
                    add_version_hashes=False,
                    project=project,
                )

                # write to disk
                io.write_dataset(
                    ds,
                    project=project,
                    output_path=new_path,
                    output_format="zarr",
                    overwrite=True,
                    compute=True,
                )


if __name__ == "__main__":
    main()
