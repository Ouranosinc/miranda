## note this needs to run using argis python 2.7
import glob
import multiprocessing as mp
import os
from os import path

import arcpy

inrep = r'H:\TRAVIS_Foret\PROJETS\PAVICS\VMshare\Projects\mddelcc_data\obs_raw\Horaire'
outrep = inrep.replace('obs_raw', 'csv_data')


def main():
    mdbs = sorted(list(glob.glob(path.join(inrep, '*.mdb'))))
    # for mdb in mdbs:
    #     export_mdb(mdb)
    n = min(len(mdbs), 3)
    pool = mp.Pool(n)

    pool.map(export_mdb, mdbs)
    pool.close()
    pool.join()


def export_mdb(mdb):
    outrep1 = path.join(outrep, path.basename(mdb).replace('.mdb', ''))
    if not path.exists(outrep1):
        os.makedirs(outrep1)
    print(mdb)
    arcpy.env.workspace = mdb
    tables = arcpy.ListTables()
    for t in tables:
        try:
            arcpy.TableToTable_conversion(t, outrep1, t + '.csv')

        except:
            print['error converting table "', t, '" to excel']


if __name__ == '__main__':
    main()
