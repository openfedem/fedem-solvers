<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# Source code organization

The FEDEM solver sources are organized into several sub-folders reflecting
the modularisation into compilation units for libraries and executables.
The important modules are as follows (listed by order of significance).

* [vpmSolver](vpmSolver) - The dynamics solver including the control system solver.
* [vpmReducer](vpmReducer) - The FE part reducer and direct linear solver.
* [vpmStress](vpmStress) - FE recovery modules
  (stress, modes, strain gage and strain coat recovery).
* [vpmCommon](vpmCommon) - Library of common Fortran modules used by all the
  three modules above. This also includes some linear algebra handling
  (equation solver encapsulation).
* [vpmUtilities](vpmUtilities) - Some low-level general utility modules.
* [Femlib](Femlib) - The Finite Element library.

## Source code documentation

We use [doxygen](https://www.doxygen.nl/) to extract documentation
from comments in the source files themselves.
See [here](https://openfedem.github.io/fedem-solvers/solver/)
for the extracted documentation (for dynamics solver and FE part reducer
and their dependencies only).
This documentation is updated automatically when a new release is tagged,
provided there are changes in the source code/comments compared with
the previous release and that these changes affects the documentation.
