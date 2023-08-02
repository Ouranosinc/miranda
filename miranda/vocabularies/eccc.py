"""Definition lists of variables from ECCC for each type of archive."""

# For more information see the ECCC Technical Documentation

__all__ = [
    "DLY",
    "DLY02",
    "DLY03",
    "DLY04",
    "DLY12",
    "DLY13",
    "DLY21",
    "DLY44",
    "HLY",
    "HLY01",
    "HLY01_RCS",
    "HLY03",
    "HLY10",
    "HLY15",
    "HLY21",
    "MLY",
    "MLY04",
]

# Hourly Data

HLY01 = []
HLY01.extend(list(range(71, 123)))  # Hourly variables
HLY01.extend([209, 210])  # Wind character and gust speed
HLY01.extend(list(range(219, 231)))  # Cloud layers
HLY01.append(244)  # Precipitation type
HLY01.append(260)  # Freezing fog

HLY01_RCS = HLY01.copy()
HLY01_RCS.extend(
    list(range(262, 281))
)  # Reference Climate Surface (RCS) weather stations

HLY03 = []
HLY03.extend(list(range(123, 133)))  # Hourly rainfall
HLY03.extend([160, 161])

HLY10 = []
HLY10.extend(list(range(61, 69)))  # Sunshine
HLY10.extend([133, 169, 170, 171, 172])  # Solar radiation

HLY15 = [69, 70, 76, 156]  # Wind

HLY21 = [123]  # Fischer/Porter precipitation

HLY = list(set(HLY01 + HLY01_RCS + HLY03 + HLY10 + HLY15 + HLY21))

# Daily Data

DLY02 = []
DLY02.extend(list(range(1, 26)))  # Daily variables
DLY02.append(157)  # Direction of extreme gust
DLY02.append(179)  # Daily bright sunshine

DLY03 = []
DLY03.extend(list(range(124, 133)))
DLY03.extend([160, 161])

DLY04 = DLY02.copy()

DLY12 = []
DLY12.extend(list(range(134, 151)))  # Soil temperatures

DLY13 = list(range(151, 156))  # Pan evaporation

DLY21 = [12]  # Precipitation
DLY21.extend(list(range(127, 133)))  # Precipitation over time
DLY21.append(161)  # Most precipitation in 25 hours

DLY44 = []
DLY44.extend([1, 2, 3])  # Temperature
DLY44.extend(list(range(10, 18)))  # Precipitation

DLY = list(set(DLY02 + DLY03 + DLY04 + DLY12 + DLY13 + DLY21 + DLY44))

# Monthly data

MLY04 = []
MLY04.extend(list(range(26, 39)))  # Days with variables
MLY04.extend(list(range(39, 61)))  # Means of variables
MLY04.append(158)  # Direction of extreme gust

MLY = list(set(MLY04))
