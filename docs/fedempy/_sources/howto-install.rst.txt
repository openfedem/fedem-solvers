Installing fedempy
==================

To install the latest released ``fedempy`` package, run the following command in a console window::

    pip install https://github.com/openfedem/fedem-solvers/releases/download/#FEDEM_TAG#/fedempy-#FEDEMPY_VERSION#.tar.gz

Before using the ``fedempy`` modules, the following environment variables need to be set:

| **FEDEM_REDUCER** = path to the Fedem reducer shared object library (`fedem_reducer_core.dll`)
| **FEDEM_SOLVER** = path to the Fedem solver shared object library (`fedem_solver_core.dll`)
| **FEDEM_MDB** = path to the Fedem mechanism database shared object library (`FedemDB.dll`)

The current version of ``fedempy`` is #FEDEMPY_VERSION# and is compatible with **Fedem #FEDEM_VERSION#**.
