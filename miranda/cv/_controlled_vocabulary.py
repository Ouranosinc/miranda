import pyessv

from ..ecmwf._era5 import ERA5_PROJECT_NAMES  # noqa

__all__ = [
    "ACTIVITIES",
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
    "WCRP_INSTITUTIONS",
]


# Controlled Vocabularies supplied via pyessv

CMIP5 = pyessv.WCRP.CMIP5  # noqa
CMIP6 = pyessv.WCRP.CMIP6  # noqa
CORDEX = pyessv.WCRP.CORDEX  # noqa
CORDEX_ADJUST = pyessv.WCRP.CORDEX_ADJUST  # noqa

# Institutes

CMIP5_INSTITUTES = [f.raw_name for f in CMIP5.institution_id.terms]
CMIP6_INSTITUTES = [f.raw_name for f in CMIP6.institution_id.terms]
CORDEX_INSTITUTES = [f.raw_name for f in CORDEX.institute.terms]
CORDEX_ADJUST_INSTITUTES = [f.raw_name for f in CORDEX_ADJUST.institute.terms]

WCRP_INSTITUTIONS = list()
WCRP_INSTITUTIONS.extend([CMIP6_INSTITUTES, CMIP5_INSTITUTES, CORDEX_INSTITUTES])

INSTITUTIONS = list()
INSTITUTIONS.extend(WCRP_INSTITUTIONS)

# Models

CMIP5_MODELS = [f.raw_name for f in CMIP5.model.terms]
CMIP6_MODELS = [f.raw_name for f in CMIP6.source_id.terms]

CORDEX_MODELS = list()
CORDEX_MODELS.extend([f.raw_name for f in CORDEX.rcm_name.terms])
CORDEX_MODELS.extend([f.raw_name for f in CORDEX.rcm_model.terms])
CORDEX_MODELS.extend(
    ["UQAM-CRCM5", "UQAM-CRCM5-SN"]
)  # Needed for internal-ish CORDEX data

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
