import logging
from argparse import ArgumentParser

import numpy as np
import pandas as pd
from xclim.core.units import convert_units_to, pint_multiply, str2pint

if __name__ == "__main__":
    argparser = ArgumentParser(description="Convert snow data to netCDF.")
    argparser.add_argument(
        "-v", "--verbose", help="Increase verbosity of the script.", action="store_true"
    )
    argparser.add_argument(
        "-o",
        "--output",
        help="Output netCDF filename.",
        default="MELCC_snow_observations.nc",
    )
    argparser.add_argument(
        "file", help="Path to the excel file with the manual snow observations."
    )
    args = argparser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    # Stations
    logging.info("Parsing stations.")
    stations = pd.read_excel(args.file, sheet_name="Stations")
    stations = stations.rename(
        columns={
            "No": "station",
            "Nom": "station_name",
            "LAT(°)": "lat",
            "LONG(°)": "lon",
            "ALT(m)": "elevation",
            "OUVERTURE": "station_opening",
            "FERMETURE": "station_closing",
        }
    )
    statds = stations.set_index("station").to_xarray()
    stations = statds.set_coords(statds.data_vars.keys()).station
    stations["station_name"] = stations["station_name"].astype(str)
    stations.lat.attrs.update(units="degree_north", standard_name="latitude")
    stations.lon.attrs.update(units="degree_east", standard_name="longitude")
    stations.elevation.attrs.update(units="m", standard_name="height")
    stations.station_opening.attrs.update(description="Date of station creation.")
    stations.station_closing.attrs.update(description="Date of station closure.")

    # Periods
    logging.info("Parsing observation periods.")
    periods = pd.read_excel(
        args.file, sheet_name="Périodes standards", names=["start", "end", "middle"]
    )
    periods = periods[["start", "end"]].to_xarray()
    periods = (
        periods.to_array()
        .rename(variable="bnds", index="num_period")
        .drop_vars("bnds")
        .rename("period_bnds")
    )
    periods.attrs.update(
        description="Bounds of the sampling periods of the MELCC. Observations are taken manually once per period. The year of these bounds should be ignored."
    )

    # Data
    logging.info("Parsing data.")
    data = pd.read_excel(
        args.file,
        sheet_name="Données",
        names=[
            "station",
            "time",
            "snd",
            "snd_flag",
            "sd",
            "sd_flag",
            "snw",
            "snw_flag",
        ],
    )
    ds = data.set_index(["station", "time"]).to_xarray()
    ds["station"] = stations.sel(station=ds.station)
    bins = periods.dt.dayofyear
    bins = np.concatenate((bins.isel(bnds=0), bins.isel(bnds=1, num_period=[-1]) + 1))
    ds["period"] = ds.time.copy(data=np.digitize(ds.time.dt.dayofyear, bins))
    ds.period.attrs.update(description="Observational period number.")

    flag_attrs = dict(
        standard_name="status_flag",
        flag_values=[0, 1, 3, 5, 7],
        flag_meanings="nodata good estimated forced trace",
        flag_meanings_fr="sansdonnée correcte estimée forcée trace",
    )
    ds.snd.attrs.update(
        standard_name="surface_snow_thickness",
        units="cm",
        long_name="Snow depth",
        long_name_fr="Épaisseur de la neige au sol",
        melcc_code="NS000F",
        melcc_description="Épaisseur de la neige mesurée (carottier)",
    )
    ds["snd"] = convert_units_to(ds.snd, "m")
    ds["snd_flag"] = ds.snd_flag.fillna(0).astype(int)
    ds.snd_flag.attrs.update(
        long_name="Quality of snow depth measurements.",
        long_name_fr="Qualité de la mesure d'épaisseur de la neige",
        **flag_attrs,
    )
    ds.snw.attrs.update(
        standard_name="surface_snow_amount",
        units="cm",
        long_name="Snow amount",
        long_name_fr="Quantité de neige au sol",
        melcc_code="NSQ000F",
        melcc_description="Équivalent en eau de la neige mesurée (carottier)",
        description="Converted from snow water-equivalent using a water density of 1000 kg/m³",
    )
    ds["snw"] = pint_multiply(ds.snw, str2pint("1000 kg m-3"), out_units="kg m^-2")
    ds["snw_flag"] = ds.snd_flag.fillna(0).astype(int)
    ds.snw_flag.attrs.update(
        long_name="Quality of snow amount measurements.",
        long_name_fr="Qualité de la mesure de quantité de neige",
        **flag_attrs,
    )
    # Density given as a percentage of water density
    ds.sd.attrs.update(
        standard_name="surface_snow_density",
        units="%",
        long_name="Snow density",
        long_name_fr="Densité de la neige au sol",
        melcc_code="NSD000F",
        melcc_description="Densité de la neige mesurée (carottier)",
    )
    ds["sd"] = pint_multiply(ds.sd, str2pint("1000 kg m-3"), out_units="kg m^-3")
    ds["sd_flag"] = ds.sd_flag.fillna(0).astype(int)
    ds.sd_flag.attrs.update(
        long_name="Quality of snow density measurements.",
        long_name_fr="Qualité de la mesure de densité de la neige",
        **flag_attrs,
    )

    # Save
    logging.info("Saving to file.")
    ds.to_netcdf(args.output)
