from os import getenv
from pathlib import Path

from miranda.eccc import convert_ahccd

if __name__ == "__main__":
    in_files = getenv("in")
    out_files = getenv("out")

    source_files = Path(in_files)
    output_path = Path(out_files)

    source_var_gens = {
        "Generation3/Homog_daily_mean_temp_v2019/": ("tas", 3),
        "Generation3/Homog_daily_max_temp_v2019/": ("tasmax", 3),
        "Generation3/Homog_daily_min_temp_v2019/": ("tasmin", 3),
        "Generation2/Adj_Daily_Total_v2017/": ("pr", 2),
        "Generation2/Adj_Daily_Snow_v2017/": ("prsn", 2),
        "Generation2/Adj_Daily_Rain_v2017/": ("prlp", 2),
    }

    for folder, (variable, generation) in source_var_gens.items():
        convert_ahccd(
            source_files.expanduser().joinpath(folder),
            output_path,
            variable,
            generation,
        )
