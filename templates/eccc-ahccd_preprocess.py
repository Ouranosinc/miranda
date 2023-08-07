from pathlib import Path

from miranda.preprocess import convert_ahccd, merge_ahccd

in_files = Path("~/Desktop/ec_data/ahccd").expanduser()
output = Path().cwd().parent / "test"
variable = "tas"

convert_ahccd(in_files, output, variable, generation=3)
merge_ahccd(output.joinpath("tas"), output.joinpath("merged"), variable, overwrite=True)
