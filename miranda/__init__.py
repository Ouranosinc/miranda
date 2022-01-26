__author__ = "Trevor James Smith"
__email__ = "smith.trevorj@ouranos.ca"
__version__ = "0.2.0-beta"


from miranda import (
    archive,
    convert,
    decode,
    deh_melcc,
    ecmwf,
    gis,
    ops,
    remove,
    scripting,
    subset,
    utils,
)

from .archiver import archive_database
from .connect import Connection
from .data import DataBase
from .decode import metadata
from .storage import FileMeta, StorageState
