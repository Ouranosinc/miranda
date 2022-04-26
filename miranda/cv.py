import pyessv

from miranda.ecmwf._era5 import ERA5_PROJECT_NAMES  # noqa

__all__ = [
    "ACTIVITIES",
    "BIAS_ADJUST_INSTITUTIONS",
    "CMIP5_INSTITUTES",
    "CMIP5_MODELS",
    "CMIP6_ACTIVITIES",
    "CMIP6_INSTITUTES",
    "CMIP6_MODELS",
    "CORDEX_INSTITUTES",
    "CORDEX_MODELS",
    "INSTITUTIONS",
    "PROJECT_MODELS",
    "REANALYSIS_ACTIVITIES",
    "WCRP_FREQUENCIES",
    "WCRP_INSTITUTIONS",
]

# Controlled Vocabularies supplied via PyESSV

CMIP5 = pyessv.WCRP.CMIP5  # noqa
CMIP6 = pyessv.WCRP.CMIP6  # noqa
CORDEX = pyessv.WCRP.CORDEX  # noqa
CORDEX_ADJUST = pyessv.WCRP.CORDEX_ADJUST  # noqa

# Institutes

CMIP5_INSTITUTES = [f.raw_name for f in CMIP5.institute.terms]
CMIP6_INSTITUTES = [f.raw_name for f in CMIP6.institution_id.terms]
CORDEX_INSTITUTES = [f.raw_name for f in CORDEX.institute.terms]
CORDEX_INSTITUTES.append("OURANOS")  # Needed for internal-ish CORDEX data
CORDEX_INSTITUTES.append("ISU")  # Needed for validation - Not presently in PyESSV
CORDEX_ADJUST_INSTITUTES = [f.raw_name for f in CORDEX_ADJUST.institute.terms]

WCRP_INSTITUTIONS = list()
WCRP_INSTITUTIONS.extend(CMIP5_INSTITUTES)
WCRP_INSTITUTIONS.extend(CMIP6_INSTITUTES)
WCRP_INSTITUTIONS.extend(CORDEX_INSTITUTES)
WCRP_INSTITUTIONS.extend(CORDEX_ADJUST_INSTITUTES)

INSTITUTIONS = list()
INSTITUTIONS.extend(WCRP_INSTITUTIONS)

BIAS_ADJUST_INSTITUTIONS = ["OURANOS", "PCIC"]

# Models

CMIP5_MODELS = [f.raw_name for f in CMIP5.model.terms]
CMIP6_MODELS = [f.raw_name for f in CMIP6.source_id.terms]

CORDEX_MODELS = list()
CORDEX_MODELS.extend([f.raw_name for f in CORDEX.rcm_name.terms])
CORDEX_MODELS.extend([f.raw_name for f in CORDEX.rcm_model.terms])
CORDEX_MODELS.extend(
    ["UQAM-CRCM5", "UQAM-CRCM5-SN"]  # Needed for internal-ish CORDEX data
)

# Time Frequencies

CMIP5_FREQUENCIES = [f.raw_name for f in CMIP5.time_frequency.terms]
CMIP6_FREQUENCIES = [f.raw_name for f in CMIP6.frequency.terms]
CORDEX_FREQUENCIES = [f.raw_name for f in CMIP5.time_frequency.terms]

WCRP_FREQUENCIES = list()
WCRP_FREQUENCIES.extend(CMIP5_FREQUENCIES)
WCRP_FREQUENCIES.extend(CMIP6_FREQUENCIES)
WCRP_FREQUENCIES.extend(CORDEX_FREQUENCIES)
WCRP_FREQUENCIES.append("sem")  # needed for some CORDEX datasets

REANALYSIS = list()
REANALYSIS.extend(ERA5_PROJECT_NAMES)
PROJECT_MODELS = dict(
    CMIP5=CMIP5_MODELS, CMIP6=CMIP6_MODELS, CORDEX=CORDEX_MODELS, REANALYSIS=REANALYSIS
)

ACTIVITIES = list()
CMIP6_ACTIVITIES = [f.raw_name for f in CMIP6.activity_id.terms]
ACTIVITIES.extend(CMIP6_ACTIVITIES)
REANALYSIS_ACTIVITIES = ["ERA"]
ACTIVITIES.extend(REANALYSIS_ACTIVITIES)
