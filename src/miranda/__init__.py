"""Python utilities for climate data collection, conversion, and management."""

###################################################################################
# Apache Software License 2.0
#
# Copyright (c) 2019-2025, Trevor James Smith
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###################################################################################

from __future__ import annotations

__author__ = """Trevor James Smith"""
__email__ = "smith.trevorj@ouranos.ca"
__version__ = "0.6.0-dev.18"


from . import (
    convert,
    cv,
    decode,
    io,
    preprocess,
    scripting,
    structure,
    units,
    utils,
    validators,
    vocabularies,
)
from .storage import FileMeta, StorageState
