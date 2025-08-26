"""Definition lists of variables from ECCC for each type of archive."""

# For more information see the ECCC Technical Documentation

__all__ = [
    "obs_groupings",
    "obs_vocabularies",
]

obs_vocabularies = dict()

# Hourly Data

obs_vocabularies["HLY01"] = []
obs_vocabularies["HLY01"].extend(list(range(71, 123)))  # Hourly variables
obs_vocabularies["HLY01"].extend([209, 210])  # Wind character and gust speed
obs_vocabularies["HLY01"].extend(list(range(219, 231)))  # Cloud layers
obs_vocabularies["HLY01"].append(244)  # Precipitation type
obs_vocabularies["HLY01"].append(260)  # Freezing fog

obs_vocabularies["HLY01_RCS"] = obs_vocabularies["HLY01"].copy()
obs_vocabularies["HLY01_RCS"].extend(list(range(262, 281)))  # Reference Climate Surface (RCS) weather stations

obs_vocabularies["HLY03"] = []
obs_vocabularies["HLY03"].extend(list(range(123, 133)))  # Hourly rainfall
obs_vocabularies["HLY03"].extend([160, 161])

obs_vocabularies["HLY10"] = []
obs_vocabularies["HLY10"].extend(list(range(61, 69)))  # Sunshine
obs_vocabularies["HLY10"].extend([133, 169, 170, 171, 172])  # Solar radiation

obs_vocabularies["HLY15"] = [69, 70, 76, 156]  # Wind

obs_vocabularies["HLY21"] = [123]  # Fischer/Porter precipitation

# Daily Data

obs_vocabularies["DLY02"] = []
obs_vocabularies["DLY02"].extend(list(range(1, 26)))  # Daily variables
obs_vocabularies["DLY02"].append(157)  # Direction of extreme gust
obs_vocabularies["DLY02"].append(179)  # Daily bright sunshine

obs_vocabularies["DLY03"] = []
obs_vocabularies["DLY03"].extend(list(range(124, 133)))
obs_vocabularies["DLY03"].extend([160, 161])

obs_vocabularies["DLY04"] = obs_vocabularies["DLY02"].copy()

obs_vocabularies["DLY12"] = []
obs_vocabularies["DLY12"].extend(list(range(134, 151)))  # Soil temperatures

obs_vocabularies["DLY13"] = list(range(151, 156))  # Pan evaporation

obs_vocabularies["DLY21"] = [12]  # Precipitation
obs_vocabularies["DLY21"].extend(list(range(127, 133)))  # Precipitation over time
obs_vocabularies["DLY21"].append(161)  # Most precipitation in 25 hours

obs_vocabularies["DLY44"] = []
obs_vocabularies["DLY44"].extend([1, 2, 3])  # Temperature
obs_vocabularies["DLY44"].extend(list(range(10, 18)))  # Precipitation

# Monthly data

obs_vocabularies["MLY04"] = []
obs_vocabularies["MLY04"].extend(list(range(26, 39)))  # Days with variables
obs_vocabularies["MLY04"].extend(list(range(39, 61)))  # Means of variables
obs_vocabularies["MLY04"].append(158)  # Direction of extreme gust

# Groupings

obs_groupings = dict()
obs_groupings["HLY"] = list(
    set(
        obs_vocabularies["HLY01"]
        + obs_vocabularies["HLY01_RCS"]
        + obs_vocabularies["HLY03"]
        + obs_vocabularies["HLY10"]
        + obs_vocabularies["HLY15"]
        + obs_vocabularies["HLY21"]
    )
)
obs_groupings["DLY"] = list(
    set(
        obs_vocabularies["DLY02"]
        + obs_vocabularies["DLY03"]
        + obs_vocabularies["DLY04"]
        + obs_vocabularies["DLY12"]
        + obs_vocabularies["DLY13"]
        + obs_vocabularies["DLY21"]
        + obs_vocabularies["DLY44"]
    )
)
obs_groupings["MLY"] = list(set(obs_vocabularies["MLY04"]))
