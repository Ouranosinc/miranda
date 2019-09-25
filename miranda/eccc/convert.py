######################################################################
# S.Biner, Ouranos, mai 2019
#
# description/commentaires
#
# programme qui decode les donnees obtenues par Trevor d'ECCC
# les donnees sont dans le reperoire data/de_trevor/Climato
# et semblent correspondre a ce qui est decrit dans
# /home/biner/projets/donnees_obs_quotid/ref/Technical_Documentation.pdf**
#
# methodologie
#
# on scan les fichiers sources annuels en cherchant une variable et on sauve
# ce qu'on trouve dans des fichiers netcdf. On applique aussi les flags
# et on fait les changements d'unites
#
# obtenu via http://climate.weather.gc.ca/index_e.html en cliquant sur 'about the data'
#######################################################################
import logging.config
from datetime import datetime as dt
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr

from miranda.scripting import LOGGING_CONFIG
from miranda.utils import get_info_var
from miranda.utils import make_local_dirs

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["convert_hourly_ec_files"]


def convert_hourly_ec_files(
    source_files: Union[str, Path],
    output_folder: Union[str, Path],
    variable_name: str,
    missing_value: int = -9999,
):
    info = get_info_var(variable_name)
    variable_code = info["code_var"]

    # on prepare l'extraction des donnees
    titre_colonnes = "code year month day code_var ".split()
    for i in range(1, 25):
        titre_colonnes.append("D{:0n}".format(i))
        titre_colonnes.append("F{:0n}".format(i))

    # preparation du repertoire de sortie
    rep_nc = Path(output_folder).joinpath("split/{}".format(nom_var))
    make_local_dirs(rep_nc)

    # boucle sur les fichiers
    list_files = Path(source_files).rglob("HLY*.gz")
    for fichier in list_files:
        logging.info("Processing file: {}.".format(fichier))

        # on converti le fichier en dataframe
        df = pd.read_fwf(
            fichier, widths=[7, 4, 2, 2, 3] + [6, 1] * 24, names=titre_colonnes
        )

        # boucle sur les eccc
        l_codes = df["code"].unique()
        for code in l_codes:
            df_code = df[df["code"] == code]

            # on arrete si la variable n'est pas presente
            if variable_code not in df_code["code_var"].unique():
                continue

            # on fait le traitement
            logging.info("Traitement de {} pour station {}".format(nom_var, code))

            # on va chercher le dataframe de la variable
            df_var = df_code[df_code["code_var"] == variable_code].copy()

            # patch pour dealer avec les bugs dans les donnees
            if (
                (fichier == str(Path(source_files).joinpath("HLYCA_123.gz")))
                and (nom_var == "rainfall")
                and (code == "7077571")
            ):
                # on corrige la donnees du 6 juin 2014 a 1h
                ind_ligne = df_var[
                    (df_var["year"] == 2014)
                    & (df_var["month"] == 6)
                    & (df_var["day"] == 14)
                ].index

                df_var.loc[ind_ligne, "D2"] = missing_value
                df_var.loc[ind_ligne, "F2"] = "M"

            # application du masque selon la valeur manquante
            df_var = df_var.replace(missing_value, np.nan)

            # decodage des valeurs et flags
            dfd = df_var.loc[:, ["D{:0n}".format(i) for i in range(1, 25)]]
            dff = df_var.loc[:, ["F{:0n}".format(i) for i in range(1, 25)]]

            # enleve flag == nan
            dff = dff.fillna("")

            # utilise flag pour masquer valeurs
            val = np.asfarray(dfd.values)
            flag = dff.values
            mask = np.isin(flag, info["flag_manquants"])
            val[mask] = np.nan

            # traitement des unites
            val = val * info["fact_mlt"] + info["fact_add"]

            # on bati le dataarray
            dates = []
            for index, row in df_var.iterrows():
                for h in range(0, 24):
                    dates.append(dt(int(row.year), int(row.month), int(row.day), h))

            ds = xr.Dataset()
            da_val = xr.DataArray(val.flatten(), coords=[dates], dims=["time"])
            da_val = da_val.rename(nom_var)
            da_val.attrs["unites"] = info["unites"]
            da_val.attrs["id"] = code
            da_val.attrs["num_element"] = variable_code

            da_flag = xr.DataArray(flag.flatten(), coords=[dates], dims=["time"])
            ds[nom_var] = da_val
            ds["flag"] = da_flag

            # sauvegarde en fichier netcdf
            an_d = ds.time.dt.year.values[0]
            an_f = ds.time.dt.year.values[-1]
            f_nc = "{c}_{v}_{ad}_{af}.nc".format(c=code, v=nom_var, ad=an_d, af=an_f)

            ds.Conventions = "CF-1.5"

            ds.title = "Environment and Climate Change Canada (ECCC) weather eccc"
            ds.history = "{}: Merged from multiple individual station files to n-dimensional array.".format(
                dt.now().strftime("%Y-%m-%d %X")
            )
            ds.version = "v{}".format(dt.now().strftime("%Y.%M"))
            ds.institution = "Environment and Climate Change Canada (ECCC)"
            ds.source = "Weather Station data <ec.services.climatiques-climate.services.ec@canada.ca>"
            ds.references = (
                "https://climate.weather.gc.ca/doc/Technical_Documentation.pdf"
            )
            ds.comment = "Acquired on demand from data specialists at ECCC Climate Services / Services Climatiques"
            ds.redistribution = "Redistribution policy unknown. For internal use only."

            ds.to_netcdf(rep_nc.joinpath(f_nc))


if __name__ == "__main__":
    # nom_var = "tas", "hourly_rainfall"
    nom_var = "precipitation_amount"
    output = "/home/tjs/Desktop/ec_data/"
    source_data = "/home/tjs/Desktop/ec_data"
    convert_hourly_ec_files(
        source_files=source_data, output_folder=output, variable_name=nom_var
    )
