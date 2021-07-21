# run using python 3.x
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar


def main():
    rootdir = Path('/exec/logan/Projets/mddelcc_data')
    # input directory containing .csv of database tables
    inrep = rootdir.joinpath('csv_data')

    #Path to stations metadata description
    freqs = {'Horaire':inrep.joinpath("STATIONs_DDs\STATIONs_DDs_Hor\STATIONs.xls"),
             'Quotidien':inrep.joinpath("STATIONs_DDs\STATIONs_DDs_Quot\STATIONs.xls")}

    # convert each .csv to a single netcdf (all stations for individual code)
    # Note multiple codes can refer to the same variable and a single station can have multiple codes for a single variable (e.g. tmax data sources by instrument w/ time overlap))
    # This script does not try to merge these common datasources - a second treatment will attempt to create into a single unified timeseries for a given variable

    # loop over hourly and daily frequencies
    for freq, infile in freqs.items():

        stat_df = pd.read_excel(infile) # read in station metadata

        # convert station start end date columns to strings
        for fld in ['Date_Fermeture', 'Date_Ouverture']:
            stat_df[fld] = stat_df[fld].dt.strftime('%Y-%m-%d')
            stat_df[fld] = stat_df[fld].replace(np.NaN, '')

        #Uniform column name look
        stat_df.columns= stat_df.columns.str.upper()
        # stat_col_dict = {'NO_Seq_Station'.upper():'station', 'NO_STATION_CLIMATO':'station_id', 'NOM_STATION':'station_name', 'Latitude'.upper():'lat',
        #        'Longitude'.upper():'lon', 'Altitude'.upper():'elev', 'Type_Poste'.upper():'station_type', 'Date_Ouverture'.upper():'fromdate',
        #        'Date_Fermeture'.upper():'todate'}

        # Get metadata on Variable codes
        dds_df = pd.read_excel(Path(infile).parent.joinpath('DDs.xls'))

        # Overall variable irregardless of equipment / station type
        #dds_df['Variable'] = dds_df['Description'].str.split("(").str[0]

        for groupe in dds_df.GROUPE.unique():

            # multiple codes for same variable (depends on sensor type etc) ... will add this info to variable metadata below
            var_codes = dds_df.loc[dds_df.GROUPE == groupe].NOM_DD

            #convert
            for var_code in var_codes:

                csv = list(inrep.joinpath(freq).rglob(f'*{var_code}.csv'))

                if len(csv) > 1:
                    raise Exception(f'expected a single .csv for {var_code} found {len(csv)}')
                if len(csv) == 1:
                    csv = csv[0]

                    # outfile = Path(csv.parent.as_posix().replace(' ', '_').lower().replace('csv_data', 'nc_data')).joinpath(
                    #     f"{variable[0:-1].replace(' ', '_').lower()}.nc")
                    outfile = Path(csv.parent.as_posix().replace(' ', '_').lower().replace('csv_data', 'nc_data')).joinpath(csv.name.replace('.csv','.nc'))

                    var_meta = dds_df[dds_df['NOM_DD'] == var_code]
                    convert_csv2nc(csv, var_meta, stat_df, outfile)

def convert_csv2nc(csv, var_meta, stat_df, outfile):

    df = pd.read_csv(csv)
    ds_list = []  # we will make a list of xr dataset (1 per station)
    for stat in df['NO_SEQ_STATION'].unique():
        print(stat, csv.name) # TODO use logging instead

        var_meta_ds = var_meta.copy()
        var_meta_ds['NO_SEQ_STATION'] = stat
        var_meta_ds.index = var_meta_ds['NO_SEQ_STATION']
        var_meta_ds = var_meta_ds.to_xarray()

        stat_meta = stat_df[stat_df['NO_SEQ_STATION'] == stat]
        stat_meta.index = stat_meta['NO_SEQ_STATION']
        stat_meta = stat_meta.drop('NO_SEQ_STATION', axis=1).to_xarray()
        if len(stat_meta['NO_SEQ_STATION'])>1:
            stat_meta = stat_meta.isel(NO_SEQ_STATION=0).expand_dims('NO_SEQ_STATION',0)


        df1 = df.loc[df['NO_SEQ_STATION'] == stat]
        date1 = pd.to_datetime(df1['DATE'], format='%d/%m/%Y %H:%M:%S')
        df1.index = date1
        df1.index.name = 'time'
        ds = df1.to_xarray().sortby('time').drop_vars(['NO_SEQ_STATION','DATE']).assign_coords(NO_SEQ_STATION=stat)
        if 'OID' in ds.data_vars:
            ds = ds.drop_vars('OID')

        ds = ds.expand_dims('NO_SEQ_STATION',0)
        ds['PRIORITE'] = var_meta_ds['PRIORITE']
        for vv in [vv for vv in var_meta_ds.data_vars if 'PRIORITE' not in vv]:
            ds['VALEUR_DONNEE'].attrs[vv]= var_meta_ds[vv][0].astype(str).values
        # add station metadata variables
        for vv in stat_meta.data_vars:
            ds[vv] = stat_meta[vv]

        ds_list.append(ds)

    if len(ds_list)>0:
        for ii, ds in enumerate(ds_list):
            # remove duplicate times (shouldn't be present but occasionally are?)
            if any(ds.get_index("time").duplicated()):
                ds_list[ii] = ds.sel(time=~ds.get_index("time").duplicated())
        with ProgressBar():
            dsOut = xr.concat(ds_list, dim='NO_SEQ_STATION').sortby('NO_SEQ_STATION','time')
            outfile.parent.mkdir(exist_ok=True,parents=True)
            dsOut.to_netcdf(outfile, format='NETCDF4_CLASSIC')

if __name__ == '__main__':
    main()
