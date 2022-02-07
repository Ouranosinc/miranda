import logging.config
from datetime import datetime as dt
from typing import Dict, List, Tuple, Union

from miranda.scripting import LOGGING_CONFIG

__all__ = ["cf_hourly_metadata", "cf_daily_metadata", "ahccd_metadata"]

logging.config.dictConfig(LOGGING_CONFIG)


def cf_hourly_metadata(variable_code: Union[int, str]) -> Dict[str, Union[int, float]]:
    """

    Parameters
    ----------
    variable_code: Union[int, str]

    Returns
    -------
    dict
    """
    ec_hourly_variables = {
        "061": {
            "nc_units": "W m-2 h-1",
            "scale_factor": 1e6 / (60 * 60),
            "add_offset": 0,
            "long_name": "RF1 Global Solar Radiation",
            "standard_name": "solar_radiation_flux",
            "nc_name": "rf1_radiation",
        },
        "071": {
            "nc_units": "m",
            "scale_factor": 30,
            "add_offset": 0,
            "long_name": "Ceiling height of lowest layer of clouds",
            "standard_name": "ceiling_cloud_height",
            "nc_name": "ceiling_hgt",
        },
        "072": {
            "nc_units": "m",
            "scale_factor": 100,
            "add_offset": 0,
            "long_name": "Visibility",
            "standard_name": "visibility_in_air",
            "nc_name": "visibility",
        },
        "073": {
            "nc_units": "Pa",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Sea Level Pressure",
            "standard_name": "air_pressure_at_mean_sea_level",
            "nc_name": "psl",
        },
        "074": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Dew Point Temperature",
            "standard_name": "dew_point_temperature",
            "nc_name": "tds",
        },
        "075": {
            "nc_units": "degree",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Wind Direction at 2 m (U2A Anemometer) (16 pts)",
            "standard_name": "wind_direction_u2a",
            "nc_name": "wind_dir_u2a_16",
        },
        "076": {
            "nc_units": "m s-1",
            "scale_factor": 0.277777778,
            "add_offset": 0,
            "long_name": "Wind Speed (U2A Anemometer)",
            "standard_name": "wind_speed_u2a",
            "nc_name": "wind_speed_u2a",
        },
        "077": {
            "nc_units": "Pa",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Station Pressure",
            "standard_name": "atmospheric_pressure",
            "nc_name": "pressure",
        },
        "078": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Dry Bulb Temperature",
            "standard_name": "dry_bulb_temperature",
            "nc_name": "tas_dry",
        },
        "079": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Wet Bulb temperature",
            "standard_name": "wet_bulb_temperature",
            "nc_name": "tas_wet",
        },
        "080": {
            "nc_units": "%",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Relative Humidity",
            "standard_name": "relative_humidity",
            "nc_name": "hur",
        },
        "081": {
            "nc_units": "%",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Total Cloud Opacity",
            "standard_name": "cloud_albedo",
            "nc_name": "clo",
        },
        "082": {
            "nc_units": "%",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Total Cloud Amount",
            "standard_name": "cloud_area_fraction",
            "nc_name": "clt",
        },
        "089": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Freezing Rain",
            "standard_name": "freezing_rain",
            "nc_name": "freeze_rain",
        },
        "094": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Ice Pellets",
            "standard_name": "ice_pellet_presence",
            "nc_name": "ice_pellets",
        },
        "107": {
            "nc_units": "%",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Lowest cloud layer opacity",
            "standard_name": "low_type_cloud_opacity_fraction",
            "nc_name": "cloud_opac",
        },
        "108": {
            "nc_units": "%",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Lowest cloud layer amount or condition",
            "standard_name": "low_type_cloud_area_fraction",
            "nc_name": "cloud_frac",
        },
        "109": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Lowest cloud layer type",
            "standard_name": "low_type_cloud_type",
            "nc_name": "low_cloud_type",
        },
        "110": {
            "nc_units": "m",
            "scale_factor": 30,
            "add_offset": 0,
            "long_name": "Lowest cloud layer height",
            "standard_name": "low_type_cloud_height",
            "nc_name": "low_cloud_hgt",
        },
        "123": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Rainfall",
            "standard_name": "rainfall_flux",
            "nc_name": "rainfall",
        },
        "133": {
            "nc_units": "s",
            "scale_factor": 3600,
            "add_offset": 0,
            "long_name": "Sunshine",
            "standard_name": "duration_of_sunshine",
            "nc_name": "sun",
        },
        "156": {
            "nc_units": "degree",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Wind Direction at 2 m (U2A Anemometer) (36 pts)",
            "standard_name": "wind_direction_u2a",
            "nc_name": "wind_dir_u2a_36",
        },
        "262": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 00-60)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation",
        },
        "263": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 00-15)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_q1",
        },
        "264": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 15-30)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_q2",
        },
        "265": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 30-45)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_q3",
        },
        "266": {
            "nc_units": "kg m-2 s-1",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Total Precipitation (minutes 45-60)",
            "standard_name": "precipitation_flux",
            "nc_name": "precipitation_q4",
        },
        "267": {
            "nc_units": "kg m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 15)",
            "standard_name": "precipitation_amount",
            "nc_name": "precipitation_weight_q1",
        },
        "268": {
            "nc_units": "kg m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 30)",
            "standard_name": "precipitation_amount",
            "nc_name": "precipitation_weight_q2",
        },
        "269": {
            "nc_units": "kg m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 45)",
            "standard_name": "precipitation_amount",
            "nc_name": "precipitation_weight_q3",
        },
        "270": {
            "nc_units": "kg m-2",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Precipitation Gauge Weight per Unit Area (at minute 60)",
            "standard_name": "precipitation_amount",
            "nc_name": "precipitation_weight_q4",
        },
        "271": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 00-15)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed_q1",
        },
        "272": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 15-30)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed_q2",
        },
        "273": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 30-45)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed_q3",
        },
        "274": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 45-60)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed_q4",
        },
        "275": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 60)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd_q4",
        },
        "276": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 15)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd_q1",
        },
        "277": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 30)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd_q2",
        },
        "278": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow Depth (at minute 45)",
            "standard_name": "surface_snow_thickness",
            "nc_name": "snd_q3",
        },
        "279": {
            "nc_units": "degree",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Wind Direction at 2 m (minutes 50-60)",
            "standard_name": "wind_direction",
            "nc_name": "wind_dir",
        },
        "280": {
            "nc_units": "m s-1",
            "scale_factor": 0.02777778,
            "add_offset": 0,
            "long_name": "Wind Speed at 2 m (minutes 50-60)",
            "standard_name": "wind_speed",
            "nc_name": "wind_speed",
        },
    }
    code = str(variable_code).zfill(3)
    if code in ["061"]:
        raise NotImplementedError()
    try:
        variable = ec_hourly_variables[code]
        variable["missing_flags"] = "M"
        variable["least_significant_digit"] = ""
    except KeyError:
        logging.error("Hourly variable `{}` not supported.".format(code))
        raise
    return variable


def cf_daily_metadata(variable_code: Union[int, str]) -> Dict[str, Union[int, float]]:
    """

    Parameters
    ----------
    variable_code: Union[int, str]

    Returns
    -------
    dict
    """
    ec_daily_variables = {
        "001": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Daily Maximum Temperature",
            "standard_name": "air_temperature_maximum",
            "nc_name": "tasmax",
        },
        "002": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Daily Minimum Temperature",
            "standard_name": "air_temperature_minimum",
            "nc_name": "tasmin",
        },
        "003": {
            "nc_units": "K",
            "scale_factor": 0.1,
            "add_offset": 273.15,
            "long_name": "Daily Mean Temperature",
            "standard_name": "air_temperature",
            "nc_name": "tas",
        },
        "010": {
            "nc_units": "mm",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Daily Total Rainfall",
            "standard_name": "liquid_precipitation_amount",
            "nc_name": "prlptot",
        },
        "011": {
            "nc_units": "cm",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Daily Total Snowfall",
            "standard_name": "solid_precipitation_amount",
            "nc_name": "prsntot",
        },
        "012": {
            "nc_units": "mm",
            "scale_factor": 0.1,
            "add_offset": 0,
            "long_name": "Daily Total Precipitation",
            "standard_name": "precipitation_amount",
            "nc_name": "prcptot",
        },
        "013": {
            "nc_units": "m",
            "scale_factor": 0.01,
            "add_offset": 0,
            "long_name": "Snow on the Ground",
            "standard_name": "surface_snow_thickness",
            "nc_name": "sndtot",
        },
        "014": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Thunderstorms",
            "standard_name": "thunderstorm_presence",
            "nc_name": "thunder",
        },
        "015": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Freezing rain or drizzle",
            "standard_name": "freeze_rain_drizzle_presence",
            "nc_name": "freezing_rain_drizzle",
        },
        "016": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Hail",
            "standard_name": "hail_presence",
            "nc_name": "hail",
        },
        "017": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Fog or Ice Fog",
            "standard_name": "fog_ice_fog_presence",
            "nc_name": "fog_ice_fog",
        },
        "018": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Smoke or Haze",
            "standard_name": "smoke_haze_presence",
            "nc_name": "smoke_haze",
        },
        "019": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Blowing Dust or Sand",
            "standard_name": "blowing_dust_sand_presence",
            "nc_name": "blowing_dust_sand",
        },
        "020": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Blowing snow",
            "standard_name": "blowing_snow_presence",
            "nc_name": "blow_snow",
        },
        "021": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Wind speed >= 28 Knots",
            "standard_name": "wind_exceeding_28_knots",
            "nc_name": "wind_gt_28kt",
        },
        "022": {
            "nc_units": "1",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "Wind speed >= 34 Knots",
            "standard_name": "wind_exceeding_34_knots",
            "nc_name": "wind_gt_34kt",
        },
        "023": {
            "nc_units": "degree",
            "scale_factor": 10,
            "add_offset": 0,
            "long_name": "Direction of extreme gust (16 pts) to December 1976",
            "standard_name": "wind_to_direction",
            "nc_name": "gust_dir",
        },
        "024": {
            "nc_units": "m s-1",
            "scale_factor": 0.2777778,
            "add_offset": 0,
            "long_name": "Speed of extreme gust",
            "standard_name": "wind_speed_of_gust",
            "nc_name": "gust_speed",
        },
        "025": {
            "nc_units": "h",
            "scale_factor": 1,
            "add_offset": 0,
            "long_name": "UTC hour of extreme gust",
            "standard_name": "hour_of_extreme_gust",
            "nc_name": "gust_hour",
        },
    }
    code = str(variable_code).zfill(3)
    try:
        variable = ec_daily_variables[code]
        variable["missing_flags"] = "M"
        variable["least_significant_digit"] = ""
    except KeyError:
        logging.error("Daily variable `{}` not supported.".format(code))
        raise
    return variable


def ahccd_metadata(
    code: str, gen: int
) -> (Dict[str, Union[str, float]], Dict, List[Tuple[int, int]], int):
    """

    Parameters
    ----------
    code: {"dx", "dn", "dm", "dt", "ds", "dr"}
    gen: {1, 2, 3}

    Returns
    -------
    Dict[str, Union[str, float]], Dict, List[Tuple[int, int]], int
    """
    generation = {1: "First", 2: "Second", 3: "Third"}.get(gen)

    ec_ahccd_attrs = dict(
        dx=dict(
            variable="tasmax",
            units="deg C",
            standard_name="air_temperature",
            long_name="Near-Surface Maximum Daily Air Temperature",
            comment=f"ECCC {generation} Generation of Adjusted and Homogenized Temperature Data",
        ),
        dn=dict(
            variable="tasmin",
            units="deg C",
            standard_name="air_temperature",
            long_name="Near-Surface Minimum Daily Air Temperature",
            comment=f"ECCC {generation} Generation of Adjusted and Homogenized Temperature Data",
        ),
        dm=dict(
            variable="tas",
            units="deg C",
            standard_name="air_temperature",
            long_name="Near-Surface Daily Mean Air Temperature",
            comment=f"ECCC {generation} Generation of Adjusted and Homogenized Temperature Data",
        ),
        dt=dict(
            variable="pr",
            units="mm day-1",
            standard_name="precipitation_flux",
            long_name="Daily Total Precipitation",
            comment=f"ECCC {generation} Generation of Adjusted and Homogenized Precipitation Data",
        ),
        ds=dict(
            variable="prsn",
            units="mm day-1",
            standard_name="snowfall_flux",
            long_name="Daily Snowfall",
            comment=f"ECCC {generation} Generation of Adjusted and Homogenized Precipitation Data",
        ),
        dr=dict(
            variable="prlp",
            units="mm day-1",
            standard_name="rainfall_flux",
            long_name="Daily Rainfall",
            comment=f"ECCC {generation} Generation of Adjusted and Homogenized Precipitation Data",
        ),
    )
    try:
        variable = ec_ahccd_attrs[code]
        variable["missing_flags"] = "M"
        if variable["variable"].startswith("tas"):
            variable["NaN_value"] = -9999.9
            column_names = [
                "No",
                "StnId",
                "Station name",
                "Prov",
                "FromYear",
                "FromMonth",
                "ToYear",
                "ToMonth",
                "%Miss",
                "Lat(deg)",
                "Long(deg)",
                "Elev(m)",
                "Joined",
                "RCS",
            ]
            column_spaces = [(0, 5), (5, 6), (6, 8), (8, 9)]
            ii = 9
            for i in range(1, 32):
                column_spaces.append((ii, ii + 7))
                ii += 7
                column_spaces.append((ii, ii + 1))
                ii += 1
            header_row = 3

        elif variable["variable"].startswith("pr"):
            variable["NaN_value"] = -9999.99
            column_names = [
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
            column_spaces = [(0, 4), (4, 5), (5, 7), (7, 8)]
            ii = 8
            for i in range(1, 32):
                column_spaces.append((ii, ii + 8))
                ii += 8
                column_spaces.append((ii, ii + 1))
                ii += 1
            header_row = 0

        else:
            raise KeyError

        column_names = {
            col.lower()
            .split("(")[0]
            .replace("%", "pct_")
            .strip()
            .replace(" ", "_"): col
            for col in list(column_names)
        }

        if gen == 3:
            _citation = (
                "Vincent, L.A., M.M. Hartwell and X.L. Wang, 2020: A Third Generation of Homogenized "
                "Temperature for Trend Analysis and Monitoring Changes in Canada’s Climate. "
                "Atmosphere-Ocean. https://doi.org/10.1080/07055900.2020.1765728"
            )
        elif gen == 2:
            _citation = (
                "Mekis, É and L.A. Vincent, 2011: An overview of the second generation adjusted daily "
                "precipitation dataset for trend analysis in Canada. Atmosphere-Ocean 49(2), "
                "163-177 doi:10.1080/07055900.2011.583910"
            )
        else:
            msg = f"Generation '{gen}' not supported."
            raise NotImplementedError(msg)

        global_attrs = dict(
            title=f"{generation} Generation of Homogenized Daily {variable['variable']} "
            "for Canada (Updated to December 2019)",
            history=f"{dt.today().strftime('%Y-%m-%d')}: Convert from original format to NetCDF",
            type="station_obs",
            institute="Environment and Climate Change Canada",
            institute_id="ECCC",
            dataset_id=f"AHCCD_gen{gen}_day_{variable['variable']}",
            frequency="day",
            licence_type="permissive",
            licence="https:/open.canada.ca/en/open-government-licence-canada",
            citation=_citation,
        )

    except KeyError as e:
        msg = f"AHCCD variable '{code}' or generation '{gen}' not supported."
        logging.error(msg)
        raise NotImplementedError(msg) from e

    return variable, column_names, column_spaces, header_row, global_attrs
