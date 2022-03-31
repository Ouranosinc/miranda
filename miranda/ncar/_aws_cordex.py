import ast
import json
import logging
import logging.config
from json.decoder import JSONDecodeError
from pathlib import Path
from typing import Dict, List, Union

import intake
import schema
import xarray as xr
from dask.diagnostics import ProgressBar

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


_allowed_args = schema.Schema(
    {
        schema.Optional("variable"): schema.Schema(
            schema.Or(
                [
                    "hurs",
                    "huss",
                    "pr",
                    "prec",
                    "ps",
                    "rsds",
                    "sfcWind",
                    "tas",
                    "tasmax",
                    "tasmin",
                    "temp",
                    "tmax",
                    "tmin",
                    "uas",
                    "vas",
                ]
            )
        ),
        schema.Optional("frequency"): "day",
        schema.Optional("scenario"): schema.Schema(
            schema.Or(["eval", "hist", "rcp45", "rcp85", "hist-rcp45", "hist-rcp85"])
        ),
        schema.Optional("grid"): schema.Or(["NAM-22i", "NAM-44i"]),
        schema.Optional("bias_correction"): schema.Or(
            ["raw", "mbcn-Daymet", "mbcn-gridMET"]
        ),
    }
)


def cordex_aws_download(
    target_folder: Union[str, Path], *, search: Dict[str, Union[str, List[str]]]
):
    schema.Schema(_allowed_args).validate(search)

    # Define the catalog description file location.
    catalog_url = (
        "https://ncar-na-cordex.s3-us-west-2.amazonaws.com/catalogs/aws-na-cordex.json"
    )

    # Interpret the "na-cordex-models" column as a list of values.
    col = intake.open_esm_datastore(
        catalog_url,
        csv_kwargs={"converters": {"na-cordex-models": ast.literal_eval}},
    )

    col_subset = col.search(**search)

    dsets = col_subset.to_dataset_dict(
        zarr_kwargs={"consolidated": True}, storage_options={"anon": True}
    )
    logging.info(f"\nDataset dictionary keys:\n {dsets.keys()}")

    dds = list()
    for key in list(dsets.keys()):
        logging.info(f"Adding {key} to the search criteria.")
        dds.append(dsets[key])

    with ProgressBar():
        for ds in dds:
            scen = ds.attrs["experiment_id"]
            for i, member in enumerate(ds.member_id):
                for var in ds.variables:
                    if var in search["variable"]:
                        var_out = var

                new_attrs = dict()
                for key, vals in ds.attrs.items():
                    try:
                        mapped = json.loads(vals)
                        if isinstance(mapped, dict):
                            new_attrs[key] = mapped.get(str(member.values), "")
                    except JSONDecodeError:
                        new_attrs[key] = vals
                    except TypeError:
                        new_attrs[key] = vals[0]

                years, datasets = zip(*ds.isel(member_id=i).groupby("time.year"))
                for d in datasets:
                    d.attrs.update(new_attrs)

                out_folder = target_folder.joinpath(f"{member.values}_{scen}")
                out_folder.mkdir(exist_ok=True)

                file_name_pattern = f"{var_out}_{member.values}_day_{scen}_{search['grid']}_{search['bias_correction']} "

                logging.info(f"Writing out files for {file_name_pattern}.")
                paths = [
                    out_folder.joinpath(f"{file_name_pattern}_{y}.nc")
                    for y in years
                    if not out_folder.joinpath(f"{file_name_pattern}_{y}.nc").exists()
                ]
                xr.save_mfdataset(datasets, paths, format="NETCDF4_CLASSIC")
