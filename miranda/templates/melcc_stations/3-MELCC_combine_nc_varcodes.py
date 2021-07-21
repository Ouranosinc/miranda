# run using python 3.x
from pathlib import Path
import pandas as pd
from dask.diagnostics import ProgressBar
from dask.distributed import Client
import xarray as xr
import numpy as np
import warnings
import tempfile
import shutil

warnings.simplefilter('ignore')

def main():

    xr.set_options(keep_attrs=True)
    rootdir = Path("/exec/logan/Projets/mddelcc_data")
    # input directory containing .csv of database tables
    inrep = rootdir.joinpath("csv_data")

    #Path to stations metadata description
    freqs = {'Horaire': inrep.joinpath("STATIONs_DDs/STATIONs_DDs_Hor/STATIONs.xls"),
             'Quotidien':inrep.joinpath("STATIONs_DDs/STATIONs_DDs_Quot/STATIONs.xls")}

    for freq, infile in freqs.items():

        #stat_df = pd.read_excel(infile)# read
        dds_df = pd.read_excel(Path(infile).parent.joinpath('DDs.xls'))

        # Overall variable irregardless of equipment / station type
        #dds_df['Variable'] = dds_df['Description'].str.split("(").str[0]
        codes = []
        # parse variable codes from dds
        # first time we can convert a sequence of 3 characters to int stop and take prefix as varcode
        for vv in dds_df['NOM_DD']:

            for i in range(0,len(vv)-2):
                #print(vv[i:i+3])
                try:
                    int(vv[i:i+3])
                    code = vv[0:i]
                    codes.append(code)
                    break
                except:
                    continue

        dds_df['variable'] = codes
        if freq == 'Quotidien':
           dds_df = dds_df.loc[dds_df.variable != 'PI']

        for vv in dds_df['variable'].unique():
            with ProgressBar():#Client(local_dir=r'/exec/logan/dask', dashboard_address=8998):
                outnc = Path('nc_data_final').joinpath(freq, f"{vv}_{freq[0]}.nc")
                outnc.parent.mkdir(parents=True, exist_ok=True)
                print(vv)
                df = dds_df.loc[dds_df.variable == vv]
                dslist = {}
                for dd in df['NOM_DD'].unique():
                    nc = list(inrep.parent.joinpath('nc_data').rglob(f"{dd}.nc"))
                    if nc:
                        if len(nc)>1:
                            raise ValueError(f'expected a single .nc for {dd} found {len(nc)} ...')

                        chunks = dict(NO_SEQ_STATION=10,time=365*15)
                        ds = xr.open_dataset(nc[0], chunks=chunks, decode_timedelta=False)
                        if len(np.unique(ds.PRIORITE)) > 1:
                            raise ValueError(f'expected a single priority value for {dd} found {len(np.unique(ds.PRIORITE))} ...')
                        dslist[np.unique(ds.PRIORITE)[0]] = ds
                        del ds

                dslist


                for ii, prior in enumerate(sorted(list(dslist.keys()))):
                    # mask statuses that are even as nan
                    print(vv, ii, prior)
                    statuts = dslist[prior].CODE_STATUT_DONNEE.to_dataframe().CODE_STATUT_DONNEE.unique()
                    if any(statuts % 2 == 0):
                        for v in dslist[prior].data_vars:
                            if 'time' in dslist[prior][v].dims:
                                dslist[prior][v] = dslist[prior][v].where(dslist[prior].CODE_STATUT_DONNEE % 2 != 0)

                    if ii == 0:
                        dsout1 = dslist[prior]
                    else:
                        #TODO overlapping days with lower priority but better status?

                        dsout1 = xr.open_dataset(outnc, chunks=chunks)
                        # align coordinates
                        ds1 = xr.Dataset()
                        dsout = xr.Dataset()
                        for v in dsout1.data_vars:
                            dsout[v], ds1[v] = xr.align(dsout1[v], dslist[prior][v], join='outer')
                        dsout.attrs = dsout1.attrs
                        dsout1 = dsout
                        del dsout

                        # get unique status values
                        statuts = ds1.CODE_STATUT_DONNEE.to_dataframe().CODE_STATUT_DONNEE.unique()
                        # loop over - highest prior
                        statuts = sorted(statuts[~np.isnan(statuts)])
                        for statut in statuts:
                            # keep dsout values where dsout is not null, otherwise replace with ds1 values for status == statut
                            for v in dsout1.data_vars:

                                if 'time' in dsout1[v].dims:
                                    dsout1[v] = dsout1[v].where(~(dsout1.VALEUR_DONNEE.isnull()), ds1[v].where(ds1.CODE_STATUT_DONNEE==statut))



                    dsout1 = convert_dtype(dsout1)

                    with tempfile.TemporaryDirectory(dir=rootdir.as_posix()) as dir1:

                        dsout1.to_netcdf(Path(dir1).joinpath('tmp.nc'), format='NETCDF4_CLASSIC')
                        dsout1.close()
                        del dsout1
                        shutil.move(Path(dir1).joinpath('tmp.nc'), outnc)



def convert_dtype(ds):
    for v in ds.data_vars:
        if ds[v].dtype=='O':
            ds[v] = ds[v].astype(str)
    return ds



if __name__ == '__main__':
    main()
