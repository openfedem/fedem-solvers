Creating/editing Fedem models
=============================

This section gives a brief introduction on how you can use the ``fedempy``
package to generate new Fedem models through python scripting.
It relies on the :mod:`modeler` module.

The modeling methods are collected in the class ``FedemModeler``,
so to access them, start your python script by::

    from fedempy.modeler import FedemModeler
    from fedempy.enums import FmDof, FmDofStat, FmLoadType, FmType, FmVar

You also need to set the environment variable **FEDEM_MDB** to point to
the shared object library of the Fedem mechanism model database,
before executing the script.
This library is named `libFedemDB.so` on Linux (`FedemDB.dll` on Windows).

Opening a new/existing model
----------------------------

To establish a Fedem model object, use::

    myModel = FedemModeler("mymodel.fmm")

If the file "mymodel.fmm" exists, this will open that model, and any subsequent
modeling operations will alter or append objects to that model.
If you want to force opening a new empty model, specify `True` as the second
parameter, viz.::

    myModel = FedemModeler("mymodel.fmm", True)

The specified model file will then be overwritten if it already exists when
saving the model.

Saving and closing current model
--------------------------------

To save the current model, use::

    if not myModel.save():
        raise Exception({"Error": "Failed to save current model"})

To give it a new name (Save As)::

    if not myModel.save("newname.fmm"):
        raise Exception({"Error": "Failed to save to newname.fmm"})

When finished, the model should be closed to release all internal memory::

    myModel.close()

Creating objects
----------------

The ``FedemModeler`` class has several methods for creating mechanism objects
of different type. They all have a signature like
``make_<object_type> ("description", [attributes])`` and return an integer
which is (with two exceptions, see below) the base Id of the generated object.
This value can be used as a handle to that object in other modeling operations.
If the object could not be created or an error occurred,
either zero or a negative value is returned.
Your modeling script should therefore always check that the returned value is
positive before continuing.

See the :class:`modeler.FedemModeler` documentation
for a full overview of the available methods.
In the following, some example statements for creating objects are presented.

FE parts
^^^^^^^^

A FE model stored in one of the supported FE data file formats can be imported
as an FE part into the current Fedem model, by using::

    p1 = myModel.make_fe_part(fe_data_file)
    print("Created a FE Part with base Id", p1)

where `fe_data_file` is the full path of the FE data file to be imported.
No other attributes of the part can be specified with this method.
See :ref:`Editing parts` below, for how to change it's properties.

Triads
^^^^^^

To create a triad object at the global location (1.0, 2.0, 3.5), use::

    t1 = myModel.make_triad("My first triad", (1, 2, 3.5))
    print("Created a Triad with base Id", t1)

If the triad should be attached to a part, the base Id of that part may be
specified using the optional `on_part` argument, as follows::

    t2 = myModel.make_triad("My second triad", (1, 2, 3.5), on_part=part_id)
    print("Created a Triad with base Id", t2, "on Part", part_id)

This assumes that the coordinates provided match a nodal point in the FE model
of the part specified.
If the triad should be attached to ground, the base Id of the Reference plane
is specified instead, which usually equals 2, e.g.::

    t3 = myModel.make_triad("My third triad", (1, 2, 3.5), on_part=2)
    print("Created a grounded Triad with base Id", t3)

Finally, if the nodal point of the triad to be created is known,
you may specify that instead of the coordinates, viz.::

    t4 = myModel.make_triad("My fourth triad", node=node_id, on_part=part_id)
    print("Created a Triad with base Id", t4, "on Part", part_id)

where `node_id` is the Id of a nodal point in the FE model of the part `part_id`.

See :ref:`Editing triads` below, for how to further modify the properties of
created triads.

Beam elements
^^^^^^^^^^^^^

To create a string of three beam elements connected to the four triads,
`t1`, `t2`, `t3` and `t4`, you can use::

    beams = myModel.make_beam("My beams", [t1, t2, t3, t4], prop_id)
    print("Created beam elements with base Ids", beams)

This method returns a list of base Id values for the created beam elements
(or `None` if an error occured). The last parameter in the above call is
the base Id of a cross section property, which is created by::

    mat_id = myModel.make_beam_material("Steel", (7850, 2.1e11, 0.3))
    prop_id = myModel.make_beam_section("Pipe", mat_id, (0.5, 0.45))

or::

    prop_id = myModel.make_beam_section("General", 0, section_data)

The first variant above creates a pipe cross section with outer diameter 0.5
and inner diameter 0.45, and connected to a material object with mass density
7850.0, Young's modulus 2.1e11 and Poisson's ratio 0.3.
The second variant creates a generic cross section, where ``section_data``
is a list of up to 10 cross section property values::

    section_data = [EA, EIy, EIz, GIt, rhoL, RhoIp, GAy, GAsz, sy, sz]

If you specify less than 10 values, the remaining values will be assumed
equal to zero.

Joints
^^^^^^

To attach a triad to ground using a rigid joint, use::

    joint1 = myModel.make_joint("Fixed", FmType.RIGID_JOINT, t1)

To connect two triads via a revolute joint::

    joint2 = myModel.make_joint("Hinge", FmType.REVOLUTE_JOINT, t2, t3)

The first triad specified (t2) will then be the dependent joint triad
and the second triad (t3) will be the independent triad. You can also specify
``FmType.BALL_JOINT`` and ``FmType.FREE_JOINT`` as joint type.

To connect triads via a cylindric joint, use::

    joint3 = myModel.make_joint ("Cylindric", FmType.CYLINDRIC_JOINT, t0, [t1, t2, ..., tn])

After the dependent triad (t0), a list of independent triads (t1, t2, ..., tn) is specified.
The first two list items (t1 and t2) are taken as the start and end position of the joint,
and the subsequent triads (t3, ...), will by the in-between triads.
The latter must lie on a straight line through the start and end triads, otherwise they won't be taken into account.

To create a prismatic joint, use the enum value ``FmType.PRISMATIC_JOINT`` instead.

See :ref:`Editing joints` below, for how to further modify the properties of
created joints.

Springs and Dampers
^^^^^^^^^^^^^^^^^^^

To create an axial spring with a piece-wise linear stiffness function,
the following will work::

    spr1 = myModel.make_spring("My first spring", (t1, t2),
                               xy=[[-0.1, 10.0], [0.0, 0.1], [0.1, 10.0]],
                               extrapol_type="FLAT")

    spr2 = myModel.make_spring("My second spring", (t3, t4), init_Stiff_Coeff=1000.0)

The first example above creates an axial spring connected to triad (base Id)
`t1` and `t2`, with a constant stiffness equal to 10.0 outside the interval
[-0.1,0.1], and a V-shaped stiffness function in between with minimum value
0.1 at zero spring deflection.
The second example creates a spring with a constant stiffness.

It is also possible to create several axial springs in one go by specifying
a list of `(int, int)` tuples, and referring an already existing spring
stiffness function using the `fn` keyword argument, as follows::

    spr = myModel.make_spring("My springs", [(t1, t2), (t3, t4), (t5, t6)], fn=spr_func)

where `spr_func` is the base Id of an existing spring stiffness function.

To create axial dampers, there exists a method `make_damper` with a similar
set of arguments as the `make_spring` method. That is, you can create one
or several dampers with a piece-wise linear, or constant, damping coefficient, e.g.::

    dmp = myModel.make_damper("My damper", (t3, t4), init_Damp_Coeff=100.0)

which creates an axial damper between the triads `t3` and `t4`
with a constant damping coefficient of 100.0.

Refer to the documentation of :meth:`modeler.FedemModeler.make_spring`
and :meth:`modeler.FedemModeler.make_damper` for an overview of all the keyword
arguments that may be used for these two methods.

External loads
^^^^^^^^^^^^^^

To create an external load on a triad, acting in the positive global
Z-direction, with a time-dependent magnitude, you can use::

    ldir = (0, 0, 1)  # Positive Z-direction
    load = myModel.make_load("Sine", FmLoadType.FORCE, t3, ldir, "1E6*sin(5*x)")
    print("Created an external load with base Id", load)

where you also can use ``FmLoadType.TORQUE`` as second parameter if a torque
load is wanted instead.
The last argument is a string with a math expression giving the load magnitude
as function of time (here represented by the variable "x"), in this case
a sinusoidal function with amplitude 1000000.0 and angular frequencey 5.0.

Alternatively, you may specify a general function as the load magnitude::

    load = myModel.make_load("My load", FmLoadType.FORCE, t3, ldir, fn=funcId)

where `funcId` is the user Id of an existing general function,
see :ref:`General functions` below.

Sensors
^^^^^^^

To create a sensor measuring the Z-displacement at a triad, use the following::

    s1 = myModel.make_sensor("Displacement", t3, FmVar.POS, FmDof.TX)
    print("Created sensor", s1)

where the third parameter can be any of ``FmVar.POS``, ``FmVar.LOCAL_VEL``,
``FmVar.GLOBAL_VEL``, ``FmVar.LOCAL_ACC``, ``FmVar.GLOBAL_ACC``,
``FmVar.LOCAL_FORCE`` and ``FmVar.GLOBAL_FORCE``, whereas the fourth parameter
can be ``FmDof.TX``, ``FmDof.TY``, ``FmDof.TZ``, ``FmDof.RX``, ``FmDof.RY``
or ``FmDof.RZ``.

To create a relative sensor between two triads, you can use::

    s1 = myModel.make_sensor("Relative displacement", (t3, t2), FmVar.POS, FmDof.TX)
    print("Created relative sensor", s2)

where the third parameter can be either ``FmVar.POS``, ``FmVar.VEL``, or
``FmVar.ACC``.

Note: This method returns the user Id of the created sensor - not its base Id.

General functions
^^^^^^^^^^^^^^^^^

To create a general function of time, any of the following can be used::

    # Polyline
    f1 = myModel.make_function("My func 1", xy=[[0,0], [1,1], [2,3], [3,0.5]], extrapol_type="FLAT")
    # Polyline from file
    f2 = myModel.make_function("My func 2", filename="data.asc", ch_name="Force")
    # Sine
    f3 = myModel.make_function("My func 3", frequency=1.23, amplitude=2.5)
    # Math expression
    f4 = myModel.make_function("My func 4", expression="1.2+3.0*x^2")
    # Constant function
    f5 = myModel.make_function("My func 5", value=3.5)
    # Linear function
    f6 = myModel.make_function("My func 6", slope=8.13)
    # Ramp function
    f7 = myModel.make_function("My func 7", start_val=0.3, start_ramp=1.46, slope=8.13)
    # Limited Ramp function
    f8 = myModel.make_function("My func 8", start_val=0.3, start_ramp=1.46, end_ramp=3.48, slope=8.13)
    # External function (no arguments)
    f9 = myModel.make_function("My func 9")

What type of function to create is determined by the presence of keywords
in the function argument list, as follows:

* `xy` : Polyline
* `filename` : Polyline from file
* `frequency` : Sine
* `expression` : Math expression
* `value` : Constant
* `slope` : Linear
* `start_ramp` : Ramp
* `end_ramp` : Limited Ramp

If no keywords (except the function name) are specified,
an external function (whose value is assigned directly through the method
:meth:`solver.FedemSolver.set_ext_func`) will be created. The other
function types may take other arguments in addition to those shown above,
refer to the method documentation :meth:`modeler.FedemModeler.make_function`
for a full overview.

Note: This method returns the user Id of the created function - not its base Id.

See :ref:`Editing functions` below, for how to further modify the properties of
created functions.

Strain rosettes
^^^^^^^^^^^^^^^

To create a strain rosette on an FE part, you can either specify the 3-4 node numbers to connect the
strain rosette element to, or the coordinates of 3 or 4 spatial points if the node numbers are not
known. It will then search for and use the closest node for each point::

    # Using FE node numbers
    r1 = myModel.make_strain_rosette("Gage A", part_id,
                                     nodes=[121, 122, 123],
                                     direction=(0, 1, 0))
    # Using spatial point coordinates
    r2 = myModel.make_strain_rosette("Gage B", part_id,
                                     pos=[(-1.702537, -0.5171, 1.702752),
                                          (-1.649538, -0.5171, 1.658139),
                                          (-1.630592, -0.5171, 1.73142)],
                                     direction=(0, 1, 0))
    print("Created strain rosettes", [r1, r2], "on FE part", part_id)

You may also specify other keywords, see :meth:`modeler.FedemModeler.make_strain_rosette`
for the full documentation of this method.

User-defined elements
^^^^^^^^^^^^^^^^^^^^^

User-defined elements can be included in a Fedem model, if you specify the path to
the plugin shared object library containing your element implementation
when creating the ``FedemModeler`` object, e.g.::

    myModel = FedemModeler("mymodel.fmm", True, "/usr/local/lib/libMyElmPlugin.so")

Please refer to the Fedem User's Guide for details on how the create a plugin library
for user-defined elements.
With this, you can create a string of three 2-noded elements connected to the four triads,
`t1`, `t2`, `t3` and `t4`, using::

    elms = myModel.make_udelm("My elements", [t1, t2, t3, t4], alpha1=0.03, alpha2=0.05)
    print("Created user-defined elements with base Ids", elms)

Currently, only two-noded elements are supported in ``fedempy``.
The path to the plugin library will be stored in the created Fedem model file.
Therefore, there is no need to specify it when :ref:`Solving Fedem models`
through a ``FmmSolver`` object.

Modifying existing objects
--------------------------

The ``FedemModeler`` class also has some methods for modifying existing objects
in the current model. They all have a signature like
``edit_<object_type> (base_id, **kwarg)`` and return the bool value True on
success, otherwise False. The `**kwarg` argument represents a varying list of
`keyword=value` pairs with the properties to assign to the object.

Editing triads
^^^^^^^^^^^^^^

To change the position of a triad with base Id `tid`, use the following::

    if not myModel.edit_triad(tid, Tx=1.0, Ty=0.2, Tz=3.4, Rx=30, Ry=10, Rz=5):
        print(" *** Failed to move Triad", tid)

The values specified are considered as offsets to the current position.
Thus, you can leave out those coordinate directions which should not change.
The rotational values (Rx, Ry, Rz) are Euler-ZYX angles (in degrees).
That is, first the Rz rotation is applied, then Ry and finally Rx.
If you need to rotate in a different order, that can be achieved by multiple
`edit_triad` calls.

To adjust the DOF status of a triad, use something like::

    if not myModel.edit_triad(tid, constraints={
                              "Tx" : FmDofStat.FIXED,
                              "Ty" : FmDofStat.PRESCRIBED,
                              "Tz" : FmDofStat.FREE_DYN,
                              "Rx" : FmDofStat.FIXED,
                              }):
        print(" *** Failed to constrain Triad", tid)

In this example, the triad is fixed in X-translation and X-axis rotation,
prescribed in Y-translation, and fixed during initial equilibrium only
in Z-translation. The last two DOFs (Ry and Rz) remain free.

If you need to constrain all DOFs in a triad, you can alternatively
use the keyword "All" to shorten the statement, viz.::

    if not myModel.edit_triad(tid, constraints={"All" : FmDofStat.FIXED}):
        print(" *** Failed to constrain Triad", tid)

This will then be equivalent to attaching the triad to ground.

To assign a constant load and/or prescribed motion to a triad, you can do::

    if not myModel.edit_triad(tid, load={"Tz" : 1000.0, "Ry" : 123.4},
                              motion={"Ty" : 0.01}):
        print(" *** Failed to assign load/motion to Triad", tid)

This will assign constant loads in the 3rd and 5th local DOF,
and a prescribed motion in the 2nd local DOF of the triad.

To assign a non-constant load or motion, just specify the user Id of the
general function defining the load magnitude instead of the constant value.
For instance, the following will assign a sinusoidal load in the `Tx`-DOF::

    lid = myModel.make_function("My load", frequency=12.5, amplitude=1000.0)
    if not myModel.edit_triad(tid, load={"Tx" : lid}):
        print(" *** Failed to assign load to Triad", tid)

The convention is that an integer value in the `load` and `motion` dictionary
argument is assumed to be the user Id of an existing general function, whereas
a real value is taken as the constant load/motion magnitude to be assigned.

Editing joints
^^^^^^^^^^^^^^

For joints, you can edit the same properties as shown for triads above,
except that the DOF status now can also be set to ``FmDofStat.SPRING`` or
``FmDofStat.SPRING_DYN``. The latter is the same as the former, except that
the DOF is kept fixed during initial equilibrium and (optionally) during
eigenvalue analysis. For joint DOFs with either of these status codes,
you can then assign constant `spring` and `damper` properties as well as
stress-free length change function, as follows::

    if not myModel.edit_joint(jid, spring={"Tx" : 1000.0, "Ry" : 1234.5},
                              damper={"Tx" : 100.0, "Ry" : 222.2},
                              length={"Tx" : len_id}):
        print(" *** Failed to assign spring/damper properties to joint", jid)

where `jid` is the base Id of the joint to modify and `len_id` is the user Id
of an existing general function defining the stress-free length change of the
`Tx`-DOF of the joint.

Editing parts
^^^^^^^^^^^^^

To change the position of a part with base Id `pid`, use the following::

    if not myModel.edit_part(pid, Tx=1.0, Ty=0.2, Tz=3.4, Rx=30, Ry=10, Rz=5):
        print(" *** Failed to move Part", pid)

The interpretation of the keywords Tx,Ty,...,Rz is here similar as for triads,
as explained above in :ref:`Editing triads`.

In addition, the structural damping and some reduction options can be changed
using the edit_part method, viz.::

    if not myModel.edit_part(pid, alpha1=0.001, alpha2=0.03,
                             component_modes=20, consistent_mass=True):
        print(" *** Failed to change properties for Part", pid)

Editing functions
^^^^^^^^^^^^^^^^^

The method `make_function` described above in :ref:`General functions`
will make a general function of time, by default.
To change the argument to other response variables, you can use the following::

    if not myModel.edit_function(fid, t1, FmVar.POS, FmDof.TX):
        print(" *** Failed to change argument of Function", fid)

Except for the first argument, which here is the user Id of the general function
to modify, this method takes the same set of arguments as the `make_sensor`
method discussed in :ref:`Sensors` above.

Using tags as object identifiers
--------------------------------

Each object in a Fedem model is assigned a unique base Id when it is created.
This is a positive integer value which is returned by the `make`-methods,
and can be used to refer to existing objects in other statements creating
or modifying objects. However, it is often more convenient to use a
user-defined `tag` to refer to an object, or a group of objects.

For this purpose, all `make`-methods accept the keyword `tag` for assigning a tag,
e.g., for triads::

    myModel.make_triad("My first Triad", (1.0, 0.0, 0.0), tag="T1")
    myModel.make_triad("My second triad", (2.0, 0.0, 0.0), tag="T2")
    myModel.make_triad("My third triad", (3.0, 1.5, 0.0), tag="T3")

Then, to change their properties, specify a string instead of the base Id::

    if not myModel.edit_triad("T.", constraints={"Tz" : FmDofStat.FIXED}):
        print(" *** Failed to constrain Triads")

The string may contain a regular expression and will expand into all objects
with a matching tag. The above example will therefore constrain the three triads
"T1", "T2" and "T3" in the Z-axis direction.

Example
-------

A sample python script using ``modeler`` to generate a simple model is provided
`here <https://github.com/openfedem/fedem-mdb/blob/main/test/fedempy/test_modeler.py>`_.
Another script that generates and solves the model which is used in the
`Car Suspension <https://github.com/openfedem/fedem-solver-tests/tree/main/TimeDomain/CarSuspension>`_ regression test is
`available here <https://github.com/openfedem/fedem-mdb/blob/main/test/fedempy/test_SLA.py>`_.

No-code modeling
================

This section gives a brief introduction on how you can use ``fedempy``
to create/edit Fedem models without the need of writing any python code yourself.
Instead, the model definition is encoded in a YAML-formatted input file.
This input file is parsed and converted into an equivalent Fedem model file (`.fmm`),
through the use of the :mod:`yaml_parser` module.

YAML input file syntax
----------------------

TODO: Descripe the file format here, listing the available keywords, etc.

Creating/editing a Fedem model with YAML input
----------------------------------------------

When you have finished the YAML input file, e.g., "myModel.yaml", execute the following command to process it::

    python -m fedempy.yaml_parser --input-file myModel.yaml --solve

The option `--solve` will execute the dynamics solver on the generated model.

Sample YAML input file
----------------------

`See here <https://github.com/openfedem/fedem-mdb/blob/main/test/fedempy/models/02-loader.yaml>`_
for a sample YAML input file, which will create the classical Loader model.
