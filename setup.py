# -*- coding: utf-8 -*-
"""Setup file for pyQD, the split-operator TDSE package."""


from setuptools import setup

VERSION = "0.1.0"
NAME = "plasmonicmeep"
LICENSE = "GPLv3"
DESCRIPTION = (
    "Script collection to calculate field enhancements of plasmonic nanostructures."
)
URL = "https://gitlab.com/theoretical-chemistry-jena/quantum-dynamics/plasmon-meep"
EMAIL = "fabian.droege@uni-jena.de"
AUTHOR = "Fabian G. Dröge"
REQUIRES_PYTHON = ">=3.7.0"
REQUIRED = [
    "numpy",
    "matplotlib",
    "h5py",
    "pandas",
    "joblib"
]
ENTRYPOINTS = {
    "console_scripts": [
        "plas-meep = plasmonicmeep.fdtd:main",
        "plas-field = plasmonicmeep.calculate_2d_field:main",
        "plas-vis = plasmonicmeep.visualise:main",
    ]
}
CLASSIFIERS = [
    # Trove classifiers
    # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: Implementation :: CPython",
]

setup(
    name=NAME,
    description=DESCRIPTION,
    author=AUTHOR,
    author_email=EMAIL,
    url=URL,
    install_requires=REQUIRED,
    python_requires=REQUIRES_PYTHON,
    packages=[NAME],
    entry_points=ENTRYPOINTS,
    classifiers=CLASSIFIERS,
    zip_safe=False,
    license=LICENSE,
    include_package_data=True,
    version=VERSION,
)
