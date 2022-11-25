#!/usr/bin/env python
import os
import re
import sys

from packaging import version
from setuptools import find_packages, setup

NAME = "miranda"
VERSION = "0.3.0"
DESCRIPTION = "Python utilities for climate data collection and management"
KEYWORDS = "climate meteorology archiving collection NetCDF"
URL = "https://github.com/Ouranosinc/miranda"
AUTHOR = "Trevor James Smith"
AUTHOR_EMAIL = "smith.trevorj@ouranos.ca"
REQUIRES_PYTHON = ">=3.7.0"
LICENSE = "Apache Software License 2.0"

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist bdist_wheel upload")
    sys.exit()

requirements = list()
with open("requirements.txt") as req:
    for dependency in req.readlines():
        if dependency == "pint" and version.parse(sys.version) < version.parse("3.8"):
            dependency = "pint<0.20"
        requirements.append(dependency)

full_requirements = list()
with open("requirements_full.txt") as dev:
    for dependency in dev.readlines():
        full_requirements.append(dependency)

dev_requirements = list()
with open("requirements_dev.txt") as dev:
    for dependency in dev.readlines():
        dev_requirements.append(dependency)
dev_requirements.extend(full_requirements)

docs_requirements = [
    "furo",
    "ipython",
    "nbsphinx",
    "pandoc",
    "sphinx",
    "sphinx_codeautolink",
    "sphinx_copybutton",
]

readme = open("README.rst").read()
doclink = """
Documentation
-------------

The full documentation is at https://miranda.readthedocs.io/en/latest."""
history = open("HISTORY.rst").read().replace(".. :changelog:", "")

hyperlink_replacements = {
    r":issue:`([0-9]+)`": r"`GH/\1 <https://github.com/Ouranosinc/miranda/issues/\1>`_",
    r":pull:`([0-9]+)`": r"`PR/\1 <https://github.com/Ouranosinc/miranda/pull/\1>`_",
    r":user:`([a-zA-Z0-9_.-]+)`": r"`@\1 <https://github.com/\1>`_",
}
for search, replacement in hyperlink_replacements.items():
    history = re.sub(search, replacement, history)

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=readme + "\n\n" + doclink + "\n\n" + history,
    long_description_content_type="text/x-rst",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    packages=find_packages(),
    package_dir={"miranda": "miranda"},
    include_package_data=True,
    install_requires=requirements,
    extras_require={
        "docs": docs_requirements,
        "full": full_requirements,
        "dev": dev_requirements,
    },
    license=LICENSE,
    zip_safe=False,
    keywords="climate meteorology archiving collection NetCDF",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
    project_urls={
        "About Ouranos": "https://www.ouranos.ca/en/",
        "Changelog": "https://miranda.readthedocs.io/en/stable/history.html",
        "Issue tracker": "https://github.com/Ouranosinc/miranda/issues",
    },
)
