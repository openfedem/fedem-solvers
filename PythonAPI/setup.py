# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Setuptools package definition
"""

from setuptools import find_packages, setup

# The build pipeline will update the version.txt file
with open("version.txt") as fd:
    VERSION = fd.read().rstrip()
# "A Better Pip Workflow": https://kennethreitz.org/essays/2016/02/25/a-better-pip-workflow
with open("requirements-to-freeze.txt") as fd:
    REQUIRES = list(fd)


setup(
    name="fedempy",
    version=VERSION,
    description="Python wrapper for the FEDEM modeler and solvers",
    author="Knut Morten Okstad",
    author_email="kmo@openfedem.org",
    url="openfedem.org",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    license="Apache 2.0",
    install_requires=REQUIRES,
    classifiers=["Intended Audience :: Internal", "Programming Language :: Python"],
    test_suite="test",
)
