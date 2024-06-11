Solving Fedem models
====================

This section gives a brief introduction on how you can use the ``fedempy``
package to execute the Fedem dynamics solver through python scripting.
It relies on the :mod:`fmm_solver` module.

Prerequisites
-------------

Before you start, you need to set the environment variable **FEDEM_SOLVER**
to point to the shared object library of the Fedem dynamics solver.
This library is named `libfedem_solver_core.so` on Linux
and `fedem_solver_core.dll` on Windows.
You also need to set the environment variable **FEDEM_MDB** to point to
the shared object library of the Fedem mechanism model database.
This library is named `libFedemDB.so` on Linux (`FedemDB.dll` on Windows).

The methods for conducting dynamics simulation of a Fedem model are collected
in the class ``FmmSolver`` which is a sub-class of ``FedemSolver``.
To access this class, start your python script by::

    from fedempy.fmm_solver import FmmSolver

Starting a new simulation
-------------------------

To start the solver on an existing Fedem model "mymodel.fmm", use::

    mySolver = FmmSolver("mymodel.fmm")

This will parse the model file, write out some key model size parameters, and
generate the necessary input files for the dynamics solver in the appropriate
directory structure. The solver will then start and halt after the setup stage,
during which the mechanism model is initialized in memory. If an error condition
occurs during this stage, an exception is raised and the script will abort.

You may also perform this in two separate operations, like::

    mySolver = FmmSolver()

    status = mySolver.start("mymodel.fmm", False, True)
    if status < 0:
        print(" *** Something bad happened, check fedem_solver.res", status)

This way, the script will not abort if the solver fails, and you can take the
proper action by checking the return value `status` instead. You may then also
reuse the FmmSolver object in subsequent simulation within the same script.

If the model contains FE parts, they will be reduced if needed, if you specify
`True` as the third parameter in the call to the ``start`` method. Otherwise,
it is assumed the FE models have already been reduced and the solver will fail
if the reduced matrix files are not found.

You also need to set the environment variable **FEDEM_REDUCER** to point to
the shared object library of the Fedem FE part reducer for this to work.
This library is named `libfedem_reducer_core.so` on Linux
and `fedem_reducer_core.dll` on Windows.

Solving
-------

There are many methods in ``FmmSolver`` via its parent class ``FedemSolver``
for running through the time stepping of the simulation.
See the :class:`fmm_solver.FmmSolver` and :class:`solver.FedemSolver`
documentation for a full overview of the available methods.

The easiest way to just step through the simulation would be as follows::

    function_id = 1  # User Id of a function measuring the response variable
    while mySolver.solve_next():
        res = mySolver.get_function(function_id)
        print("Time = ", mySolver.get_current_time(), " response =", res)
    if mySolver.solver_done() == 0 and mySolver.ierr.value == 0:
        print("   * Time step loop finished, solver is closed")
    else:
        print(" *** Solver failed, see fedem_solver.res", mySolver.ierr.value)

The while loop will continue as long there are more time steps defined,
and no error condition occurred. In the above, the current simulation time
and a response variable extracted using ``get_function`` is printed during
the simulation. Here you can insert any other processing instead.

If you only need to run through the model without accessing any response
variables, the above code block can be replaced by one line::

     mySolver.solve_all("mymodel.fmm", True, True)

This will also reduce the FE parts, if needed, before the solver is started,
see :meth:`fmm_solver.FmmSolver.solve_all`.

Saving/closing current model
----------------------------

To close model after the simulation has finished, use::

    mySolver.close_model(True)

This will first update the model file with references to the latest simulation
results in the model database, and then release the model from core memory.
You may also close the model without saving it (e.g., if the simulation failed)
by specifying `False` as parameter.

If you used the ``solve_all`` method, you don't need to invoke ``close_model``
aftwards as this is included in the former method.

Solving models with external functions
--------------------------------------

In digital twin applications, the input loads and/or motions in a Fedem model
are defined by means of External Functions. These are function objects that are
assigned their value by an external process, for instance by using the method
:meth:`solver.FedemSolver.set_ext_func` within the time integration loop.

During the development of a digital twin model, it is often convenient to take
the external function values from a file instead, such that the model can be
solved directly from the Fedem Desktop, or by using the ``solve_all`` method.
This can be done by specifying the following additional option for the
dynamics solver::

    -externalfuncfile <datafile.asc>

where `<datafile.asc>` is the name of a multi-column ASCII file containing the
external function values as columns. The first column is assumed to contain the
time of each step and is not used. The columns need to be labeled with the `tag`
of the external function objects, by specifying the following as a
comment line before the numerical data::

    #DESCRIPTION <func1-tag> <func2-tag> <func3-tag> ...

where each `<func#-tag>` identifies the respective column.
The order of the columns in the file is arbitrary, and it may also contain some
columns that are not used, since a search will be performed for each function.

For instance, assume you created a model using :class:`modeler.FedemModeler`
and included the following external functions::

    myModel.make_function("Input function A", tag="FuncA");
    myModel.make_function("Input function B", tag="FuncB");

Then a `datafile.asc` containing the following will work::

    #DESCRIPTION Func_A unused Func_B
    0.0          1.2345  0.123 2.3456
    0.1          2.2345  0.234 3.3456
    ...

Note: It is assumed that the time steps of the simulation match those of
the specified file, and such that a new line of data is read and the
function values are updated prior to each time step of the simulation.
No interpolation is performed if the times do not match. If the file contains
fewer data lines than the total number of time steps in the simulation,
all external functions will remain constant and equal to the values of
the last line for the remaining steps throughout the simulation.
