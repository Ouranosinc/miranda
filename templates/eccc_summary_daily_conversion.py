from pathlib import Path

from miranda.eccc import daily_summaries_to_netcdf, extract_daily_summaries

if __name__ == "__main__":
    d = extract_daily_summaries(Path().cwd().joinpath("testdata"), file_suffix="*csv*")
    daily_summaries_to_netcdf(d, Path().cwd().joinpath("processed"))
