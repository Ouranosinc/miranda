import json
import logging.config
import os
import re
from pathlib import Path
from typing import Optional, Tuple, Union
import datetime
import pandas as pd
import xarray as xr
from miranda.convert import rechunk_reanalysis

from miranda.scripting import LOGGING_CONFIG
from miranda.units import u
from _utils import variable_conversion, load_json_data_mappings
logging.config.dictConfig(LOGGING_CONFIG)
home = os.path.expanduser('~')

def main(project=None, input_folder=None, output_folder=None, output_format=None):

    var_attrs = load_json_data_mappings(project=project)['variable_entry']

    for vv in var_attrs.keys():
        year = datetime.date.today().year
        for year in range(1980, year):
            for month in range(1,13):
                ncfiles = sorted(list(input_folder.glob(f"{year}{str(month).zfill(2)}*.nc")))
                outfile_root = ncfiles[0].stem.split()
                ds_all = None
                for nc in ncfiles:
                    print(nc.name)
                    ds = xr.open_dataset(nc, chunks='auto')
                    ds1 = xr.Dataset(attrs=ds.attrs)
                    ds1[vv] = ds[vv]

                    if ds_all is None:
                        ds_all = ds1
                    else:
                        ds_all = xr.concat([ds_all, ds1], dim='time')
                del ds_all.attrs['License']
                del ds_all.attrs['product']
                ds_all.attrs['Remarks'] = ds_all.attrs['Remarks'].replace('Variable names', 'Original variable names')
                ds_all = variable_conversion(ds_all, project=project, output_format='zarr')
                if output_format == "netcdf":
                    output_folder.mkdir(exist_ok=True)
                    out = output_folder / f"{ncfiles[0].stem}.nc"
                elif output_format == "zarr":
                    output_path = output_folder / "temp"
                    output_path.mkdir(exist_ok=True, parents=True)
                    out = output_path / f"{outfile_root}.zarr"
                #_write_data(ds_all, outfile)
            #rechunk_reanalysis(project=project, ds=ds_all, variables=variables, output_folder=outrep, target_chunks=dict(time=-1, lon=50,lat=50))

def _write_data(ds=None, outfile=None):
    ds
    outfile

if __name__ == '__main__':
    project = "rdrs-v2.1"
    kwargs = dict(project = project,
                  input_folder = Path(home).joinpath('RDRS_v2.1'),
                  output_folder = Path(home).joinpath('tmpout',project),
                  output_format = "zarr"
                  )


    main(**kwargs)
