"""
Copyright 2019-2022 Trevor James Smith and Ouranos Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

__author__ = "Trevor James Smith"
__email__ = "smith.trevorj@ouranos.ca"
__version__ = "0.3.0"


from . import (
    archive,
    convert,
    cv,
    decode,
    eccc,
    ecmwf,
    ncar,
    remote,
    scripting,
    structure,
    units,
    utils,
    validators,
)
from .data import DataBase
from .storage import FileMeta, StorageState
