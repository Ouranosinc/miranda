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
# etape 1)
#
#    on scan les fichiers sources annuels en cherchant une variable et on sauve
#    ce qu'on trouve dans des fichiers netcdf. On applique aussi les flags
#    et on fait les changements d'unites
#
# etape 2)
#
#    on rassemble les fichiers netcdf des differentes stations en un
#    seul fichier netCDF.
#
# Notes:
#
# En traitant les fichiers d'entree a l'étape 1, j'ai réalisé que dans
# certains fichiers il y a plusieurs entrées (2 au max je pense) pour
# une même variable/journée/station. J'ai donc ajouté l'option de
# garder la première ou la deuxième entrée (variable keep_double =
# 'first' ou 'last'). Une comparaison avec les données obtenues par
# Bruno suggère que l'option keep_double='start' est probablement la
# meilleur.
#
#
# **obtenu via http://climate.weather.gc.ca/index_e.html en cliquant sur
#   'about the data')
from .ec_aggregate_stations import *
from .ec_convert_stations import *
