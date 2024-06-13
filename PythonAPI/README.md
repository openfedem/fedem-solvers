<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# The FEDEM python API

This folder contains the sources of the python wrapper for FEDEM, the `fedempy` package.
The source code resides in the sub-folder [src/fedempy](src/fedempy) which is
compiled into the python package `fedempy` by the [setup.py](setup.py) script.
In addition, some test drivers are placed in the sub-folder [PythonAPITests](PythonAPITests)
which are invoked as regression/integration tests by the cmake-based build system
for the `fedem-solvers` project.

See [here](https://openfedem.github.io/fedem-solvers/fedempy/)
for the extracted source code documentation of the python code,
which is generated using the [Sphinx](https://www.sphinx-doc.org) tool.
That page also contains installation and end-user documentation for `fedempy`.
This documentation is updated automatically when a new release is tagged,
provided there are changes in the source code/comments compared with
the previous release and that these changes affects the documentation.

## Local build of python package

Run the command

    python setup.py sdist

This will create the archive `fedempy-<VERSION>.tar.gz`
where `<VERSION>` is the content of the file [version.txt](version.txt).
To install the built package, use the command:

    pip install dist/fedempy-<VERSION>.tar.gz

For doing development on the `fedempy` package, you may choose to install
from the sources in editable mode instead, i.e.,

    pip install --editable .

## Local build of source code documentation

Run the bash script [doc/make.sh](doc/make.sh) to rebuild the html-documentation
using the `sphinx-build` tool. It will generate the html-files in the subfolder
`doc/build/html`.

## Static code analysis

You can perform some static code checks locally by running the following commands:

    python -m isort --check --diff setup.py src test
    python -m black --config=pyproject.toml --check --diff setup.py src test
    python -m pylint setup.py src test
    python -m mypy --config-file mypy.ini setup.py src test
