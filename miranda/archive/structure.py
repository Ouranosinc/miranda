import datetime
import logging
import os
import shutil

import netCDF4

logging.basicConfig(
    filename="{}_cordex-na.log".format(
        datetime.datetime.strftime(datetime.datetime.now(), "%Y-%m-%d")
    ),
    level=logging.INFO,
)

path_origin = "/expl7/smith/CORDEX/CORDEX-NA"
path_move = "/expl7/smith/CORDEX/CORDEX-NA/structured"

rcm_institutions = {
    "CanRCM4": "CCCMA",
    "CRCM5": "UQAM",
    "CRCM5-UQAM": "UQAM",
    "HIRHAM5": "DMI",
    "RCA4": "SMHI",
    "RegCM4": "ISU",
    "WRF": "NCAR",
}

driving_institutions = {
    "CanESM2": "CCCMA",
    "CNRM-CM5": "CNRM-CERFACS",
    "CNRM-CM5-2": "CNRM-CERACS",
    "EC-EARTH": "ICHEC",
    "ERAINT": "ECMWF",
    "GEMatm-Can": "UQAM",
    "GEMatm-MPI": "UQAM",
    "GFDL-ESM2M": "NOAA-GFDL",
    "HadGEM2-ES": "MOHC",
    "MPI-ESM-LR": "MPI-M",
    "MPI-ESM-MR": "MPI-M",
}

for root, dirs, files in os.walk(path_origin):
    for nc_file in files:
        decode_name = nc_file.split(".")
        if decode_name[-1] != "nc":
            continue
        var_name = decode_name[0]
        experiment = decode_name[1]
        if experiment == "eval":
            experiment = "evaluation"
        elif experiment == "hist":
            experiment = "historical"
        driving_data = decode_name[2]
        if driving_data == "ERA-Int":
            driving_data = "ERAINT"
        rcm = decode_name[3]
        time_frequency = decode_name[4]
        domain = decode_name[5]
        bias_correction = decode_name[6]
        experiment_tag = "{}_{}-{}_{}".format(
            domain, driving_institutions[driving_data], driving_data, experiment
        )
        nc = netCDF4.Dataset(os.path.join(root, nc_file), "r")
        try:
            ensemble_member = nc.driving_model_ensemble_member
        except AttributeError:
            ensemble_member = "r1i1p1"
            logging.info("Guessing r1i1p1 for {}".format(os.path.join(root, nc_file)))

        nc.close()
        move_path = os.path.join(
            path_move,
            "CORDEX-NA",
            rcm_institutions[rcm],
            rcm,
            experiment_tag,
            time_frequency,
            "atmos",
            ensemble_member,
            bias_correction,
            var_name,
        )
        if not os.path.isdir(move_path):
            os.makedirs(move_path)
        shutil.move(os.path.join(root, nc_file), os.path.join(move_path, nc_file))
