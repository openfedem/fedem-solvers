# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
This Python test basically does the same as `test_extfunc.py`
but now using the window dts operator as the simulation driver.

The test uses two environment variables:
    FEDEM_SOLVER = Full path to the solver shared object library
    TEST_DIR = Full path to parent folder of the solver tests in the build tree
"""

from os import environ
from math import sin, cos
from pandas import DataFrame, concat, read_csv
from fedempy.dts_operators import window
from test_utils import compare_lists


def F1(x):
    """
    First input function.
    """
    return 1.0e7 * sin(x)


def F2(x):
    """
    Second input function.
    """
    return 5.0e3 * (cos(2.5 * x) - 1.0)


class DTScontext:
    """
    Dummy DTS container class.
    """
    def __init__(self):
        self._window = None
        self._state = None

    @property
    def window(self):
        return self._window

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, value):
        self._state = value


# Define the state container
dts = DTScontext()

# Start the fedem solver in the $TEST_DIR/Cantilever-extFunc folder.
wrkdir = environ["TEST_DIR"] + "/Cantilever-extFunc"
print("\n#### Running dynamics solver in", wrkdir)

# Read reference values into a DataFrame
refdata = read_csv(wrkdir + "/TipPosition.asc", sep="\t", skiprows=6)

# Create the input data
nstep0 = 100
nsteps = nstep0 + nstep0
inputs = { "F1" : [], "F2" : [] }
tindex = []
for i in range(1,1+nsteps):
    t = refdata.iloc[i, 0]
    tindex.append(t)
    inputs["F1"].append(F1(t))
    inputs["F2"].append(F2(t))
df = DataFrame(inputs,index=tindex)

# Set up keyword dictionary for the simulation driver
slvopt = { "fco": "Setup.fco", "fsifile": "Model.fsi" }
kwargs = { "lib_dir": wrkdir, "solver_options": slvopt }
kwargs.update({ "output_ids": [3, 4, 5], "use_state": True })

# Run the simulation over the first nstep0 steps
outputs1 = window.run(df.iloc[:nstep0, :].copy(), dts, **kwargs)
print(f"\nHere are the outputs for the first {nstep0} steps:\n", outputs1)

# Run the simulation over next nstep0 steps
outputs2 = window.run(df.iloc[nstep0:, :], dts, **kwargs)
print(f"\nHere are the outputs for the next {nstep0} steps:\n", outputs2)

# Concatenate the two output frames
outputs = concat([outputs1, outputs2])
outputs.insert(0, "Time", refdata.iloc[1:1+nsteps, 0].values)
print(f"\nHere are all {nsteps} steps:\n", outputs)

# Check against reference values
ctol = 1.0e-8
ierr = 0
for i in range(nsteps):
    t = refdata.iloc[i+1, 0]
    ierr += compare_lists(t, outputs.iloc[i, :].values, refdata.iloc[i+1, :].values, ctol)
if ierr > 0:
    print(f" *** Detected {ierr} discrepancies, test not passed")
    exit(ierr)
else:
    print(f"   * All response values match, checked {nsteps} steps")
