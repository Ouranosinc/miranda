from pathlib import Path

import pandas as pd

from miranda.eccc import convert_ahccd_fwf_files

testsrc = dict()
meta_src = dict()
testsrc["tasmax"] = Path(
    r"H:\TRAVIS_Foret\PROJETS\5-Data\ECCC_CRD\AHCCD_daily\Homog_daily_max_temp_v2019"
)
meta_src["tasmax"] = list(
    Path(r"H:\TRAVIS_Foret\PROJETS\5-Data\ECCC_CRD\AHCCD_daily").glob(
        "Homog_Temperature_Stations.xls"
    )
)
testsrc["tasmin"] = Path(
    r"H:\TRAVIS_Foret\PROJETS\5-Data\ECCC_CRD\AHCCD_daily\Homog_daily_min_temp_v2019"
)
meta_src["tasmin"] = meta_src["tasmax"]
testsrc["tas"] = Path(
    r"H:\TRAVIS_Foret\PROJETS\5-Data\ECCC_CRD\AHCCD_daily\Homog_daily_mean_temp_v2019"
)
meta_src["tas"] = meta_src["tasmax"]
testsrc["pr"] = Path(
    r"H:\TRAVIS_Foret\PROJETS\5-Data\ECCC_CRD\AHCCD_daily\Adj_Daily_Total_v2017"
)
meta_src["pr"] = list(
    Path(r"H:\TRAVIS_Foret\PROJETS\5-Data\ECCC_CRD\AHCCD_daily").glob(
        "Adj_Precipitation_Stations.xls"
    )
)
testsrc["prsn"] = Path(
    r"H:\TRAVIS_Foret\PROJETS\5-Data\ECCC_CRD\AHCCD_daily\Adj_Daily_Snow_v2017"
)
meta_src["prsn"] = meta_src["pr"]
testsrc["prlp"] = Path(
    r"H:\TRAVIS_Foret\PROJETS\5-Data\ECCC_CRD\AHCCD_daily\Adj_Daily_Rain_v2017"
)
meta_src["prlp"] = meta_src["pr"]
vari_code = dict(tasmax="dx", tasmin="dn", tas="dm", pr="dt", prsn="ds", prlp="dr")
vari_attrs = dict(
    tasmax=dict(
        units="deg C",
        standard_name="air_temperature",
        long_name="Near-Surface Maximum Daily Air Temperature",
        comment="ECCC Second Generation of Homogenized Temperature Data",
    ),
    tasmin=dict(
        units="deg C",
        standard_name="air_temperature",
        long_name="Near-Surface Minimum Daily Air Temperature",
        comment="ECCC Second Generation of Homogenized Temperature Data",
    ),
    tas=dict(
        units="deg C",
        standard_name="air_temperature",
        long_name="Near-Surface Daily Mean Air Temperature",
        comment="ECCC Second Generation of Homogenized Temperature Data",
    ),
    pr=dict(
        units="mm day-1",
        standard_name="precipitation_flux",
        long_name="Total Precipitation",
        comment="ECCC Second Generation of Homogenized Precipiation Data",
    ),
    prsn=dict(
        units="mm day-1",
        standard_name="snowfall_flux",
        long_name="Snowfall",
        comment="ECCC Second Generation of Homogenized Precipiation Data",
    ),
    prlp=dict(
        units="mm day-1",
        standard_name="rainfall_flux",
        long_name="Rainfall",
        comment="ECCC Second Generation of Homogenized Precipiation Data",
    ),
)
outpath = Path(
    r"H:\TRAVIS_Foret\PROJETS\PAVICS\VMshare\CCCS-Portal\CRD_data\AHCCD_daily\netcdf\separate_stations"
)
for variable in testsrc.keys():
    outpath.joinpath(variable).mkdir(parents=True, exist_ok=True)
    if "tas" in variable:
        metadata = pd.read_excel(meta_src[variable][0], header=3)
        metadata.columns = [
            "Prov",
            "Station name",
            "stnid",
            "beg yr",
            "beg mon",
            "end yr",
            "end mon",
            "lat (deg)",
            "long (deg)",
            "elev (m)",
            "stns joined",
            "RCS",
            "%Miss",
        ]
        cols_specs = [(0, 5), (5, 6), (6, 8), (8, 9)]
        ii = 9
        for i in range(1, 32):
            cols_specs.append((ii, ii + 7))
            ii += 7
            cols_specs.append((ii, ii + 1))
            ii += 1
    else:
        metadata = pd.read_excel(meta_src[variable][0], header=3)
        metadata.columns = [
            "Prov",
            "Station name",
            "stnid",
            "beg yr",
            "beg mon",
            "end yr",
            "end mon",
            "lat (deg)",
            "long (deg)",
            "elev (m)",
            "stns joined",
        ]
        cols_specs = [(0, 4), (4, 5), (5, 7), (7, 8)]
        ii = 8
        for i in range(1, 32):
            cols_specs.append((ii, ii + 8))
            ii += 8
            cols_specs.append((ii, ii + 1))
            ii += 1
        for index, row in metadata.iterrows():
            if isinstance(row["stnid"], str):
                metadata.loc[index, "stnid"] = metadata.loc[index, "stnid"].strip(" ")
    for ff in testsrc[variable].glob("*d*.txt"):
        print(ff.name)
        stid = ff.name.replace(vari_code[variable], "").split(".txt")[0]

        try:
            metadata_st = metadata[metadata["stnid"] == int(stid)]
        except ValueError:
            metadata_st = metadata[metadata["stnid"] == stid]

        if len(metadata_st) == 1:
            dsOut = convert_ahccd_fwf_files(
                ff,
                metadata_st,
                variable,
                cols_specs=cols_specs,
                attrs=vari_attrs[variable],
            )
        else:
            raise Exception("metadata not found")

        outfile = outpath.joinpath(variable, ff.name.replace(".txt", ".nc"))
        dsOut.to_netcdf(outfile)
