# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Python wrapper for the native fedem modeler.

This module provides functionality for creating new Fedem models,
and for modifying existing ones.
"""

from ctypes import POINTER, byref, c_bool, c_char_p, c_double, c_int
from os import environ, path

from fedempy.enums import FmDof, FmDofStat, FmType, FmVar
from fedempy.fmm import FedemModel, _convert_char, _convert_int, _convert_real


def _convert_bool(arg):
    """
    Converts a boolean from Python to C type.
    """
    if arg is None:
        return c_bool(False)

    if isinstance(arg, bool):
        return c_bool(arg)

    if isinstance(arg, c_bool):
        return arg

    raise TypeError(f"Expected {bool} or {c_bool}, got {type(arg)}.")


def _convert_int_array(arg):
    """
    Converts an array of integers from Python to C type.
    """
    if arg is None:
        return c_int(0), None

    if isinstance(arg, (list, tuple)):
        argc = len(arg)
        if argc > 0:
            argv = (c_int * argc)()
            argv[:] = arg
        else:
            argv = None
        return c_int(argc), argv

    if hasattr(arg, "_length_") and getattr(arg, "_type_", None) is c_int:
        return c_int(len(arg)), arg

    raise TypeError(f"Expected {list} or {c_int}, got {type(arg)}.")


def _convert_float_array(arg):
    """
    Converts an array of floats from Python to C type.
    """
    if arg is None:
        return c_int(0), None

    if isinstance(arg, (list, tuple)):
        argc = len(arg)
        if argc > 0:
            argv = (c_double * argc)()
            argv[:] = arg
        else:
            argv = None
        return c_int(argc), argv

    if hasattr(arg, "_length_") and getattr(arg, "_type_", None) is c_double:
        return c_int(len(arg)), arg

    raise TypeError(f"Expected {list} or {c_double}, got {type(arg)}.")


def _extract_dofs(arg):
    """
    Extracts 6 DOF values into a c_double array from a dictionary.
    """
    if not isinstance(arg, dict):
        return None

    has_dofs = False
    dof_tags = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]
    dof_vals = (c_double * 6)(*([0.0] * 6))
    for key, value in arg.items():
        if key in dof_tags:
            dof_vals[dof_tags.index(key)] = _convert_real(value)
            has_dofs = True

    if has_dofs:
        return dof_vals

    return None  # Return nothing if no DOF tags were found


def _extract_polyline(arg):
    """
    Extracts data for a polyline function from a dictionary.
    """
    coord = arg["xy"]
    n_ = c_int(len(coord))
    x_ = (c_double * len(coord))()
    y_ = (c_double * len(coord))()
    x_[:] = [c_double(c[0]) for c in coord]
    y_[:] = [c_double(c[1]) for c in coord]
    extrap = ["NONE", "FLAT", "LINEAR"]
    e_type = arg.get("extrapol_type", extrap[0])
    if e_type in extrap:
        e_ = c_int(extrap.index(e_type))
    else:
        e_ = c_int(0)

    return n_, x_, y_, e_


class FedemModeler(FedemModel):
    """
    This subclass of :class:`fmm.FedemModel` adds some basic modeling methods.

    Parameters
    ----------
    model_file : str, default=None
        Absolute path to the fedem model file to open
    force_new : bool, default=False
        If True, any existing model will be overwritten on save
    plugins : str or list of str, default=None
        Plugin libraries with user-defined elements and functions
    name : str, default=None
        Description of the new model, not used if opening existing model

    Methods
    -------
    open:
        Opens the specified model file
    save:
        Saves current model into the specified model file
    close:
        Closes the currently open model
    make_triad:
        Creates a triad at specified location or node
    make_beam:
        Creates a string of beam elements
    make_beam_section:
        Creates a beam cross section property object
    make_beam_material:
        Creates a material property object
    make_spring:
        Creates a spring element
    make_damper:
        Creates a damper element
    make_joint:
        Creates a joint object
    make_load:
        Creates an external load object
    make_function:
        Creates a general function object
    make_sensor:
        Creates a sensor object
    make_fe_part:
        Creates an FE part
    make_generic_part:
        Creates a generic part
    make_strain_rosette:
        Creates a strain rosette on an FE part
    make_udelm:
        Creates a string of user-defined elements
    edit_triad:
        Modifies an existing triad
    edit_part:
        Modifies an existing FE part
    edit_joint:
        Modifies an existing joint
    edit_function:
        Modifies an existing function
    solver_setup:
        Setting dynamics solver parameters
    """

    def __init__(self, model_file=None, force_new=False, plugins=None, name=None):
        """
        Constructor.
        Optionally opens a Fedem model, if a model file name is specified.
        """
        fedem_lib = environ.get("FEDEM_MDB")
        if not fedem_lib:
            print("\n *** Environment variable FEDEM_MDB is not set!", flush=True)
            return  # This will probably cause crash later, but not here

        if isinstance(plugins, str):
            super().__init__(fedem_lib, plugins)
        elif isinstance(plugins, list) and len(plugins) > 0:
            if len(plugins) > 1:
                super().__init__(fedem_lib, plugins[0], plugins[1])
            else:
                super().__init__(fedem_lib, plugins[0])
        else:
            super().__init__(fedem_lib)

        self._fmlib.FmCreateTriad.restype = c_int
        self._fmlib.FmCreateTriad.argtypes = [
            c_char_p,
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
        ]

        self._fmlib.FmTriadOnNode.restype = c_int
        self._fmlib.FmTriadOnNode.argtypes = [c_char_p, c_int, c_int]

        self._fmlib.FmCreateBeam.restype = c_int
        self._fmlib.FmCreateBeam.argtypes = [c_char_p, c_int, c_int]

        self._fmlib.FmCreateBeamProperty.restype = c_int
        self._fmlib.FmCreateBeamProperty.argtypes = [
            c_char_p,
            c_int,
            c_int,
            POINTER(c_double),
        ]

        self._fmlib.FmCreateMaterialProperty.restype = c_int
        self._fmlib.FmCreateMaterialProperty.argtypes = [
            c_char_p,
            c_int,
            POINTER(c_double),
        ]

        self._fmlib.FmCreateSpring.restype = c_int
        self._fmlib.FmCreateSpring.argtypes = [
            c_char_p,
            c_int,
            c_int,
            c_double,
            c_bool,
            c_double,
            POINTER(c_int),
            c_int,
            POINTER(c_double),
            POINTER(c_double),
            c_int,
            c_int,
        ]

        self._fmlib.FmCreateDamper.restype = c_int
        self._fmlib.FmCreateDamper.argtypes = [
            c_char_p,
            c_int,
            c_int,
            c_bool,
            c_double,
            POINTER(c_int),
            c_int,
            POINTER(c_double),
            POINTER(c_double),
            c_int,
        ]

        self._fmlib.FmCreateJoint.restype = c_int
        self._fmlib.FmCreateJoint.argtypes = [
            c_char_p,
            c_int,
            c_int,
            POINTER(c_int),
            c_int,
        ]

        self._fmlib.FmCreateLoad.restype = c_int
        self._fmlib.FmCreateLoad.argtypes = [
            c_char_p,
            c_int,
            c_int,
            c_double,
            c_double,
            c_double,
            c_char_p,
            c_int,
        ]

        self._fmlib.FmCreateMathExprFunc.restype = c_int
        self._fmlib.FmCreateMathExprFunc.argtypes = [
            c_char_p,
            c_char_p,
            c_char_p,
            c_bool,
        ]

        self._fmlib.FmCreateExternalFunc.restype = c_int
        self._fmlib.FmCreateExternalFunc.argtypes = [c_char_p, c_char_p, c_bool]

        self._fmlib.FmCreateSineFunc.restype = c_int
        self._fmlib.FmCreateSineFunc.argtypes = [
            c_char_p,
            c_char_p,
            POINTER(c_double),
            c_bool,
        ]

        self._fmlib.FmCreateLinearFunc.restype = c_int
        self._fmlib.FmCreateLinearFunc.argstypes = [
            c_char_p,
            c_char_p,
            POINTER(c_double),
            c_bool,
        ]

        self._fmlib.FmCreatePolyFunc.restype = c_int
        self._fmlib.FmCreatePolyFunc.argtypes = [
            c_char_p,
            c_char_p,
            c_int,
            POINTER(c_double),
            POINTER(c_double),
            c_int,
            c_bool,
        ]

        self._fmlib.FmCreateDeviceFunc.restype = c_int
        self._fmlib.FmCreateDeviceFunc.argtypes = [
            c_char_p,
            c_char_p,
            c_char_p,
            c_char_p,
            c_double,
            c_bool,
            c_double,
            c_bool,
        ]

        self._fmlib.FmCreateSensor.restype = c_int
        self._fmlib.FmCreateSensor.argtypes = [
            c_char_p,
            c_char_p,
            c_int,
            c_int,
            c_int,
            c_int,
        ]

        self._fmlib.FmSetFunctionArg.restype = c_bool
        self._fmlib.FmSetFunctionArg.argtypes = [c_int, c_int, c_int, c_int, c_int]

        self._fmlib.FmLoadPart.restype = c_int
        self._fmlib.FmLoadPart.argtypes = [c_char_p, c_char_p]

        self._fmlib.FmCreatePart.restype = c_int
        self._fmlib.FmCreatePart.argtypes = [c_char_p, c_int, POINTER(c_int)]

        self._fmlib.FmCreateStrainRosette.restype = c_int
        self._fmlib.FmCreateStrainRosette.argtypes = [
            c_char_p,
            c_int,
            c_int,
            POINTER(c_int),
            POINTER(c_double),
            c_double,
            c_bool,
        ]

        self._fmlib.FmCreateUDE2.restype = c_int
        self._fmlib.FmCreateUDE2.argtypes = [c_char_p, c_int, c_int]

        self._fmlib.FmCreateAssembly.restype = c_int
        self._fmlib.FmCreateAssembly.argtypes = [c_char_p, c_int, POINTER(c_int)]

        self._fmlib.FmMoveObject.restype = c_bool
        self._fmlib.FmMoveObject.argtypes = [c_int, POINTER(c_double), c_int, c_int]

        self._fmlib.FmConstrainObject.restype = c_bool
        self._fmlib.FmConstrainObject.argtypes = [c_int, c_int, c_int]

        self._fmlib.FmAddMass.restype = c_bool
        self._fmlib.FmAddMass.argtypes = [c_int, c_int, POINTER(c_double), c_int]

        self._fmlib.FmDofProperty.restype = c_bool
        self._fmlib.FmDofProperty.argtypes = [c_int, c_int, c_int, c_double, c_int]

        self._fmlib.FmStructDamp.restype = c_bool
        self._fmlib.FmStructDamp.argtypes = [c_int, c_double, c_double]

        self._fmlib.FmReduceOpts.restype = c_bool
        self._fmlib.FmReduceOpts.argtypes = [c_int, c_int, c_bool]

        self._fmlib.FmRecoverOpts.restype = c_bool
        self._fmlib.FmRecoverOpts.argtypes = [c_int, c_int, c_bool]

        if model_file is None or not path.isfile(model_file) or force_new:
            self.fm_new(model_file, name)
        elif not self.fm_open(model_file):
            print(" *** Failed to open model file", model_file)

    def _convert_ids(self, obj_ids, obj_type):
        """
        Convenience method, to reduce cognitive complexity.
        """
        if isinstance(obj_ids, list):
            return obj_ids

        return self._convert_id(obj_ids, obj_type, True)

    def open(self, model_file):
        """
        Opens the specified model file and prints out some key model parameters.

        Parameters
        ----------
        model_file : str
            Absolute path of the Fedem model file to open

        Returns
        -------
        bool
            True on success, otherwise False
        """
        if not self.fm_open(model_file):
            print(" *** Failed to open model file", model_file)
            return False

        print("   * Model file", model_file, "successfully opened")
        print("     Number of Triads:", self.fm_count(FmType.TRIAD))
        print("     Number of Beams:", self.fm_count(FmType.BEAM))
        print("     Number of Parts:", self.fm_count(FmType.FEPART))
        print("     Total number of mechanism objects:", self.fm_count())

        return True

    def save(self, model_file=None):
        """
        Saves current model into the specified model file.

        If no model_file is given, that last opened model file is overwritten.

        Parameters
        ----------
        model_file : str, default=None
            Absolute path of the Fedem model file to save to

        Returns
        -------
        bool
            True on success, otherwise False
        """
        status = self.fm_save(model_file)
        if not status:
            print(" *** Failed to save model file", model_file)
        return status

    def close(self, save=False, remove_singletons=False):
        """
        Closes the currently open model.

        Parameters
        ----------
        save : bool, default=False
            If True, the model file is updated with the current model
        remove_singletons : bool, default=False
            If True, heap-allocated singelton objects are also released

        Returns
        -------
        bool
            True on success, otherwise False
        """
        status = self.fm_save() if save else True
        if not status:
            print(" *** Failed to save current fedem model")
            self.fm_close(True)
        else:
            self.fm_close(remove_singletons)
        return status

    def make_triad(self, name, pos=None, rot=None, node=0, on_part=0, tag=None):
        """
        Creates a new triad at specified location or nodal point.

        Parameters
        ----------
        name : str
            Description of the new triad
        pos : (float, float, float)
            Global XYZ-coordinates of new triad
        rot : (float, float, float), default=None
            Global Euler angles giving the orientation of new triad
        node : int, default=0
            FE node number to associate the triad with.
            Used only if `on_part` is specified.
        on_part : int or str, default=0
            Base Id or tag of the part that this triad should be attached to.
            You can also specify the Reference plane here,
            to create a new triad that is attached to ground.
        tag : str, default=None
            Tag to associate the created triad with

        Returns
        -------
        int
            Base Id of new triad, zero or negative on error
        """
        part_ = self._convert_id(on_part)
        if node > 0 and part_.value > 2:
            triad = self._fmlib.FmTriadOnNode(
                _convert_char(name), _convert_int(node), part_
            )

        elif len(pos) < 3:
            print(" *** Invalid Triad position", pos)
            triad = -1

        elif rot is None:
            if len(pos) > 5:
                triad = self._fmlib.FmCreateTriad(
                    _convert_char(name),
                    _convert_real(pos[0]),
                    _convert_real(pos[1]),
                    _convert_real(pos[2]),
                    _convert_real(pos[3]),
                    _convert_real(pos[4]),
                    _convert_real(pos[5]),
                    part_,
                )
            else:
                triad = self._fmlib.FmCreateTriad(
                    _convert_char(name),
                    _convert_real(pos[0]),
                    _convert_real(pos[1]),
                    _convert_real(pos[2]),
                    c_double(0),
                    c_double(0),
                    c_double(0),
                    part_,
                )

        elif len(rot) < 3:
            print(" *** Invalid Triad rotation", rot)
            triad = -1

        else:
            triad = self._fmlib.FmCreateTriad(
                _convert_char(name),
                _convert_real(pos[0]),
                _convert_real(pos[1]),
                _convert_real(pos[2]),
                _convert_real(rot[0]),
                _convert_real(rot[1]),
                _convert_real(rot[2]),
                part_,
            )

        if tag is not None and triad > 0:
            self.fm_tag_object(triad, tag)

        return triad

    def make_beam(self, name, triads, bprop=None, tag=None):
        """
        Creates a string of beam elements.

        Parameters
        ----------
        name : str
            Description of the new beam(s)
        triads : list of int or list of str
            List of base Ids or tags of the connected triads
        bprop : int or str, default=None
            Base Id or tag of beam property to use
        tag : str, default=None
            Tag to associate the created beam(s) with

        Returns
        -------
        list of int
            Base Ids of the created beams, None if error
        """
        if len(triads) < 2:
            print(" *** make_beam: At least two triads must be specified")
            return None

        base_ids = []
        for i in range(len(triads) - 1):
            base_id = self._fmlib.FmCreateBeam(
                _convert_char(name),
                self._convert_id(triads[i], FmType.TRIAD),
                self._convert_id(triads[i + 1], FmType.TRIAD),
                self._convert_id(bprop, FmType.BEAM_PROP),
            )
            if base_id < 1:
                return None
            base_ids.append(base_id)

        if tag is not None:
            self.fm_tag_object(base_ids, tag)

        return base_ids

    def make_beam_section(self, name, mat, bprops, tag=None):
        """
        Creates a beam cross section property object.

        Parameters
        ----------
        name : str
            Description of the new beam property
        mat : int or str
            Cross section type flag.
            If zero, a Generic cross section is defined.
            If str or a non-zero int, a Pipe cross section is defined,
            and the value gives the tag or base Id of the material to use.
        bprop : list of float
            List of property data
        tag : str, default=None
            Tag to associate the created beam property with

        Returns
        -------
        int
            Base Id of beam property object, zero or negative on error
        """
        nprop_, prop_ = _convert_float_array(bprops)
        bprop = self._fmlib.FmCreateBeamProperty(
            _convert_char(name),
            self._convert_id(mat, FmType.MAT_PROP),
            nprop_,
            prop_,
        )

        if tag is not None and bprop > 0:
            self.fm_tag_object(bprop, tag)

        return bprop

    def make_beam_material(self, name, mprops, tag=None):
        """
        Creates a material property object.

        Parameters
        ----------
        name : str
            Description of the new material property
        mprop : list of float
            List of property data
        tag : str, default=None
            Tag to associate the created material property with

        Returns
        -------
        int
            Base Id of material property object, zero or negative on error
        """
        nprop_, prop_ = _convert_float_array(mprops)
        mprop = self._fmlib.FmCreateMaterialProperty(_convert_char(name), nprop_, prop_)

        if tag is not None and mprop > 0:
            self.fm_tag_object(mprop, tag)

        return mprop

    def make_spring(self, name, triads, **kwargs):
        """
        Creates axial spring elements.

        Parameters
        ----------
        name : str
            Description of the new spring(s)
        triads : (str, str) or (int, int) or list of (int, int)
            Tags or base Ids of the connected triads.
            If a list of tuples is specified,
            one spring is created for each tuple.
        tag : str, default=None
            Tag to associate the created beam(s) with
        kwargs : dict
            Keyword arguments defining the spring properties.
            The following keywords are currently supported:

            * `tag` : Tag to associate the created spring(s) with
            * `constDefl`: Constant stress free deflection
            * `constLength`: Constant stress free length
            * `length`: User Id of general function defining stress free length
            * `init_Stiff_Coeff`: Constant spring stiffness coefficient
            * `fn`: Base Id of spring stiffness function
            * `xy`: List of XY-pairs giving a piece-wise linear spring stiffness
            * `extrapol_type`: String, either "NONE" (default), "FLAT" or "LINEAR"
            * `spring_characteristics`: String, either "SPR_TRA_STIFF" (default)
              or "SPR_TRA_FORCE"

        Returns
        -------
        int or list of int
            Base Id(s) of the created spring(s), None if an error occurs
        """
        n_ = c_int(0)
        x_ = None
        y_ = None
        e_ = c_int(0)

        use_constant_defl_ = c_bool(True)  # default is constant deflection
        if "constDefl" in kwargs:
            const_length_defl_ = _convert_real(kwargs["constDefl"])
        elif "constLength" in kwargs:
            const_length_defl_ = _convert_real(kwargs["constLength"])
            use_constant_defl_ = c_bool(False)
        else:
            const_length_defl_ = c_double(0)

        if "fn" in kwargs:  # Base Id of existing spring function
            sp_ = _convert_int(kwargs["fn"])
        else:
            sp_ = c_int(0)

        if "xy" in kwargs:
            # We have data of a new spring stiffness function
            n_, x_, y_, e_ = _extract_polyline(kwargs)

            charac = [
                "SPR_TRA_STIFF",  # stiffness - translational deflection
                "SPR_TRA_FORCE",  # force - translational deflection
            ]
            sp_charac = kwargs.get("spring_characteristics", charac[0])
            if sp_charac in charac:
                # Notice negative value to indicate spring function type
                sp_ = c_int(-charac.index(sp_charac))

        one_spring = isinstance(triads, tuple)
        if one_spring:
            triads = [triads]  # Only one spring is created

        base_ids = []
        for triad in triads:
            s_id = self._fmlib.FmCreateSpring(
                _convert_char(name),
                self._convert_id(triad[0], FmType.TRIAD),
                self._convert_id(triad[1], FmType.TRIAD),
                const_length_defl_,
                use_constant_defl_,
                _convert_real(kwargs.get("init_Stiff_Coeff", 0.0)),
                byref(sp_),
                n_,
                x_,
                y_,
                e_,
                _convert_int(kwargs.get("length", 0)),
            )
            if s_id > 0:
                base_ids.append(s_id)

        if len(base_ids) < len(triads):
            return None  # Failure creating at least one spring

        if "tag" in kwargs:
            self.fm_tag_object(base_ids, kwargs["tag"])

        if one_spring:
            return base_ids[0]

        return base_ids

    def make_damper(self, name, triads, **kwargs):
        """
        Creates axial damper elements.

        Parameters
        ----------
        name : str
            Description of the new damper(s)
        triads : (str, str) or (int, int) or list of (int, int)
            Tags or base Ids of the connected triads.
            If a list of tuples is specified,
            one damper is created for each tuple.
        kwargs : dict
            Keyword arguments defining the damper properties.
            The following keywords are currently supported:

            * `tag` : Tag to associate the created damper(s) with
            * `def_vel_damper`: If True, use deformational velocity
            * `init_Damp_Coeff`: Constant damping coefficient
            * `fn`: Base Id of damping coefficient function
            * `xy`: List of XY-pairs giving a piece-wise linear damping coefficient
            * `extrapol_type`: String, either "NONE" (default), "FLAT" or "LINEAR"
            * `damp_characteristics`: String, either "DA_TRA_COEFF" (default)
              or "DA_TRA_FORCE"

        Returns
        -------
        int or list of int
            Base Id(s) of the created damper(s), None if an error occurs
        """
        n_ = c_int(0)
        x_ = None
        y_ = None
        e_ = c_int(0)

        if "fn" in kwargs:  # Base Id of existing damper function
            da_ = _convert_int(kwargs["fn"])
        else:
            da_ = c_int(0)

        if "xy" in kwargs:
            # We have data of a new damping coefficient function
            n_, x_, y_, e_ = _extract_polyline(kwargs)

            charac = [
                "DA_TRA_COEFF",  # damping coefficient - translational velocity
                "DA_TRA_FORCE",  # force - translational velocity
            ]
            da_charac = kwargs.get("damp_characteristics", charac[0])
            if da_charac in charac:
                # Notice negative value to indicate damper function type
                da_ = c_int(-charac.index(da_charac))

        one_damper = isinstance(triads, tuple)
        if one_damper:
            triads = [triads]  # Only one damper is created

        base_ids = []
        for triad in triads:
            d_id = self._fmlib.FmCreateDamper(
                _convert_char(name),
                self._convert_id(triad[0], FmType.TRIAD),
                self._convert_id(triad[1], FmType.TRIAD),
                _convert_bool(kwargs.get("def_vel_damper", False)),
                _convert_real(kwargs.get("init_Damp_Coeff", 0.0)),
                byref(da_),
                n_,
                x_,
                y_,
                e_,
            )
            if d_id > 0:
                base_ids.append(d_id)

        if len(base_ids) < len(triads):
            return None  # Failure creating at least one damper

        if "tag" in kwargs:
            self.fm_tag_object(base_ids, kwargs["tag"])

        if one_damper:
            return base_ids[0]

        return base_ids

    def make_joint(self, name, joint_type, follower, followed=None, tag=None):
        """
        Creates a joint object.

        Parameters
        ----------
        name : str
            Description of the new joint
        joint_type : FmType
            Type of joint
        follower : int or str
            Base Id or tag of the dependent joint triad
        followed : int or str or list of int, default=None
            Base Id or tag of the independent joint triad(s).
            If None, the joint is connected to ground and the independent
            triad is created at the same location as the dependent triad.
            For point-to-path joints, the first two triads specified
            are taken as the end points of the glider.
        tag : str, default=None
            Tag to associate the created joint with

        Returns
        -------
        int
            Base Id of joint object, zero or negative on error
        """
        id_ = self._convert_id(follower, FmType.TRIAD)

        if isinstance(followed, list):
            nids_, ids_ = _convert_int_array(followed)
        else:
            nids_ = c_int(1)
            ids_ = (c_int * 1)(self._convert_id(followed, FmType.TRIAD))

        joint = self._fmlib.FmCreateJoint(
            _convert_char(name), _convert_int(joint_type), id_, ids_, nids_
        )

        if tag is not None and joint > 0:
            self.fm_tag_object(joint, tag)

        return joint

    def make_load(
        self, name, load_type, triad, load_dir, magnitude=None, fn=0, tag=None
    ):
        """
        Creates an external load object.

        Parameters
        ----------
        name : str
            Description of the new load
        load_type : FmLoadType
            Type of load
        triad : int or str
            Base Id or tag of triad where the load attacks
        load_dir : (float, float, float)
            Load direction vector
        magnitude : str, default=None
            Load magnitude expression
        fn : int, default=0
            User Id of load magnitude function
        tag : str, default=None
            Tag to associate the created load with

        Returns
        -------
        int
            Base Id of load object, zero or negative on error
        """
        if len(load_dir) < 3:
            print(" *** Invalid load direction vector", load_dir)
            return -1

        load = self._fmlib.FmCreateLoad(
            _convert_char(name),
            _convert_int(load_type),
            self._convert_id(triad, FmType.TRIAD),
            _convert_real(load_dir[0]),
            _convert_real(load_dir[1]),
            _convert_real(load_dir[2]),
            _convert_char(magnitude),
            _convert_int(fn),
        )

        if tag is not None and load > 0:
            self.fm_tag_object(load, tag)

        return load

    def make_function(self, name, **kwargs):
        """
        Creates a general function of time.

        The type of function to be created is determined
        by which keyword arguments are provided.
        The keyword determining the function type is below marked
        by an asterix (:sup:`*`) in each case.

        Parameters
        ----------
        name : str
            Description of the new function

        kwargs : dict
            Keyword arguments depending on function type.
            The following **function types** and `keywords`
            are currently supported:

            1. **Polyline**

              | `xy`:sup:`*`: List of `XY`-pairs giving a
                              piece-wise linear curve
              | `extrapol_type`: String, either "NONE" (default),
                                 "FLAT" or "LINEAR"

            2. **Polyline-from-file**

              | `filename`:sup:`*`: Name of file containing `XY`-pairs
              | `ch_name`: String identifying the column to use
                           for multi-column files
              | `sc_factor`: Scaling factor, default=1.0
              | `z_adjust`: If True, the `y`-values are shifted
                            such that the first value is zero, default=False
              | `v_shift`: Additional shift of the `y`-values, default=0.0

            3. **Sine**

              | `frequency`:sup:`*`: Angular frequency
              | `amplitude`: Scaling factor, default=1.0
              | `delay`: Phase shift, default=0.0
              | `mean_value`: Constant shift, default=0.0
              | `end`: Default=0.0, if greater than zero, the function
                       value is constant for `x` greater than `end`

              The Sine function therefore evaluates to::

                  f(x) = amplitude*sin(frequency*x-delay) + mean_value

              for `x`-values less than `end`, whereas ``f(x) = f(end)``
              for `x` greater than `end`.

            4. **Constant**

              | `value`:sup:`*`: The constant function value, i.e.,
                                 ``f(x) = value`` for any `x`

            5. **Linear**

              | `slope`:sup:`*`: The scaling factor, i.e., ``f(x) = slope*x``

            6. **Ramp**

              | `start_ramp`:sup:`*`: Start point of sloped function domain
              | `start_val`: Function value before `start_ramp`, default=0.0
              | `slope`: The scaling factor, default=1.0

              The Ramp function therefore evaluates to::

                  f(x) = start_val + slope*(x-start_ramp)

              | for `x`-values greater than `start_ramp`, whereas
              | ``f(x) = start_val`` for `x` less than or equal to `start_ramp`.

            7. **Limited Ramp**

              | `start_ramp`: Start point of sloped function domain, default=0.0
              | `end_ramp`:sup:`*`: End point of sloped function domain
              | `start_val`: Function value before `start_ramp`, default=0.0
              | `slope`: The scaling factor, default=1.0

              The Limited Ramp function therefore evaluates to::

                  f(x) = start_val + slope*(x-start_ramp)

              | for `x`-values in the range [`start_ramp`, `end_ramp`]
              | whereas ``f(x) = start_val`` for `x` less than `start_ramp`,
              | and ``f(x) = f(end_ramp)`` for `x` greater than `end_ramp`.

            8. **Math expression**

              | `expression`:sup:`*`: String containing the expression to use

            9. **External function**

              | No keywords

           In addition, the `tag` keyword is accepted for all function types,
           to associate a specified tag with the created function, and the
           `base_id` keyword is used to indicate if the base Id of the created
           function should be returned, instead of the default user Id.

        Returns
        -------
        int
            User Id or base Id of created function object.
            Will be 0 or negative if an error occurs.
        """

        name_ = _convert_char(name)
        tag_ = _convert_char(kwargs.get("tag", None))
        bid_ = _convert_bool(kwargs.get("base_id", False))

        if "xy" in kwargs:
            # Polyline function, get abscissa and ordinate values
            n_, x_, y_, e_ = _extract_polyline(kwargs)
            return self._fmlib.FmCreatePolyFunc(name_, tag_, n_, x_, y_, e_, bid_)

        if "frequency" in kwargs:
            # Sine function, get function parameters
            sine_par = (c_double * 5)()
            sine_par[0] = c_double(kwargs.get("frequency", 1.0))
            sine_par[1] = c_double(kwargs.get("delay", 0.0))
            sine_par[2] = c_double(kwargs.get("amplitude", 1.0))
            sine_par[3] = c_double(kwargs.get("mean_value", 0.0))
            sine_par[4] = c_double(kwargs.get("end", 0.0))
            return self._fmlib.FmCreateSineFunc(name_, tag_, sine_par, bid_)

        if "filename" in kwargs:
            # Polyline-from-file function
            return self._fmlib.FmCreateDeviceFunc(
                name_,
                tag_,
                _convert_char(kwargs.get("filename", None)),
                _convert_char(kwargs.get("ch_name", None)),
                _convert_real(kwargs.get("sc_factor", 1.0)),
                _convert_bool(kwargs.get("z_adjust", False)),
                _convert_real(kwargs.get("v_shift", 0.0)),
                bid_,
            )

        if "expression" in kwargs:
            # Math expression function
            return self._fmlib.FmCreateMathExprFunc(
                name_, tag_, _convert_char(kwargs["expression"]), bid_
            )

        if "slope" in kwargs or "value" in kwargs:
            params = (c_double * 4)(*([0.0] * 4))
            if "end_ramp" in kwargs:  # Limited ramp function
                params[0] = c_double(kwargs.get("slope", 1.0))
                params[1] = c_double(kwargs.get("start_val", 0.0))
                params[2] = c_double(kwargs.get("start_ramp", 0.0))
                params[3] = c_double(kwargs.get("end_ramp", 5.0))
            elif "start_ramp" in kwargs:  # Ramp function
                params[0] = c_double(kwargs.get("slope", 1.0))
                params[1] = c_double(kwargs.get("start_val", 0.0))
                params[2] = c_double(kwargs.get("start_ramp", 0.0))
                params[3] = params[2]
            elif "slope" in kwargs:  # Linear function
                params[0] = c_double(kwargs.get("slope", 1.0))
            else:  # Constant function
                params[1] = c_double(kwargs.get("value", 0.0))
            return self._fmlib.FmCreateLinearFunc(name_, tag_, params, bid_)

        # By default (if kwargs is empty), create an external function
        return self._fmlib.FmCreateExternalFunc(name_, tag_, bid_)

    def make_sensor(self, name, obj, var, dof=None, tag=None):
        """
        Creates a sensor object.

        Parameters
        ----------
        name : str
            Description of the new sensor
        obj : int or str or (int, int) or (str, str)
            Base Id or tag of object(s) to measure
        var : FmVar
            The variable to measure
        dof : FmDof, default=None
            Local DOF to measure if `var` is a multi-DOF quantity
        tag : str, default=None
            Tag to associate the created sensor with

        Returns
        -------
        int
            User Id of sensor object, zero or negative on error
        """
        if isinstance(obj, tuple):
            obj1_ = self._convert_id(obj[0])
            obj2_ = self._convert_id(obj[1])
        else:
            obj1_ = self._convert_id(obj)
            obj2_ = c_int(0)

        return self._fmlib.FmCreateSensor(
            _convert_char(name),
            _convert_char(tag),
            _convert_int(var),
            _convert_int(dof),
            obj1_,
            obj2_,
        )

    def make_fe_part(self, file, name=None, tag=None):
        """
        Creates an FE part.

        Parameters
        ----------
        file : str
            Absolute path to the FE data file
        name : str, default=None
            Description of the new part, use file basename if not specified
        tag : str, default=None
            Tag to associate the created part with

        Returns
        -------
        int
            Base Id of the FE part object, zero or negative on error
        """
        part = self._fmlib.FmLoadPart(_convert_char(file), _convert_char(name))

        if tag is not None and part > 0:
            self.fm_tag_object(part, tag)

        return part

    def make_generic_part(self, name, triads, **kwargs):
        """
        Creates a generic part.

        Parameters
        ----------
        name : str
            Description of the new part
        triads : list of int or list of str
            List of base Ids or tags of the connected triads
        kwargs : dict
            Keyword arguments defining some element properties.
            Currently, the following keywords are recognized:

            * `alpha1` : Mass-proportional damping coefficient
            * `alpha2` : Stiffness-proportional damping coefficient
            * `tag` : Tag to associate the created element with

        Returns
        -------
        int
            Base Id of the created generic part, zero of negative on error
        """
        n_triads = len(triads)
        if n_triads < 1:
            print(" *** make_generic_part: No triads specified")
            return 0

        n_ = c_int(n_triads)
        t_ = (c_int * n_triads)()
        tr_[:] = [self._convert_id(t, FmType.TRIAD) for t in triads]
        base_id = self._fmlib.FmCreatePart(_convert_char(name), n_, tr_)
        if base_id < 1:
            return base_id

        if "alpha1" in kwargs or "alpha2" in kwargs:
            self._fmlib.FmStructDamp(
                base_id,
                _convert_real(kwargs.get("alpha1", 0.0)),
                _convert_real(kwargs.get("alpha2", 0.0)),
            )

        if "tag" in kwargs:
            self.fm_tag_object(base_id, _convert_char(kwargs["tag"]))

        return base_id

    def make_strain_rosette(self, name, part_id, **kwargs):
        """
        Creates a strain rosette on an FE part.

        Parameters
        ----------
        name : str
            Description of the new strain rosette
        part_id : int or str
            Base Id or tag of FE part on which the rosette will be created
        kwargs : dict
            Keyword arguments defining the strain rosette properties.
            The following keywords are recognized:

            * `tag` : Tag to associate the created strain rosette with
            * `nodes`: List of nodes to connect the strain rosette to
            * `pos`: List of connection point coordinates, on the format
              `[(x1,y1,z1), (x2,y2,z2), (x3,y3,z3), (x4,y4,z4)]`.
              The last tuple is skipped if a 3-noded rosette is desired.
            * `direction`: Tuple `(x,y,z)` defining the local `X`-axis
              direction vector of the strain rosette
            * `angle` : Angle between strain gage direction and local X-axis
            * `start_at_zero` : If True, the start strains are set to zero

        Returns
        -------
        int
            Base Id of the strain rosette, zero or negative on error
        """
        name_ = _convert_char(name)
        part_ = self._convert_id(part_id, FmType.FEPART)

        nodes_ = (c_int * 4)()
        node_id = kwargs.get("nodes", None)
        node_pos = kwargs.get("pos", None)
        if isinstance(node_pos, list):
            # Nodal positions are provided, find the closest FE nodes
            nnod_ = c_int(len(node_pos))
            for idx, pos in enumerate(node_pos):
                nodes_[idx] = self.fm_get_node(part_, pos)
        elif isinstance(node_id, list):
            # The node Id's are provided directly
            nnod_ = c_int(len(node_id))
            nodes_[0 : len(node_id)] = node_id
        else:
            print(" *** Neither nodes nor pos is specified")
            return -1

        dirvec_ = (c_double * 3)()
        direction = kwargs.get("direction", None)
        if isinstance(direction, tuple):
            dirvec_[:] = direction
        else:
            dirvec_[:] = [1.0, 0.0, 0.0]

        ros = self._fmlib.FmCreateStrainRosette(
            name_,
            part_,
            nnod_,
            nodes_,
            dirvec_,
            _convert_real(kwargs.get("angle", 0.0)),
            _convert_bool(kwargs.get("start_at_zero", True)),
        )

        if ros > 0:  # Activate strain gage recovery during solve for this part
            self._fmlib.FmRecoverOpts(part_, c_int(2), c_bool(True))

        if "tag" in kwargs and ros > 0:
            self.fm_tag_object(ros, _convert_char(kwargs["tag"]))

        return ros

    def make_udelm(self, name, triads, **kwargs):
        """
        Creates a string of 2-noded user-defined elements.

        Parameters
        ----------
        name : str
            Description of the new element(s)
        triads : list of int or list of str
            List of base Ids or tags of the connected triads
        kwargs : dict
            Keyword arguments defining some element properties.
            Currently, the following keywords are recognized:

            * `alpha1` : Mass-proportional damping coefficient
            * `alpha2` : Stiffness-proportional damping coefficient
            * `tag` : Tag to associate the created element with

        Returns
        -------
        list of int
            Base Ids of the created elements, None if error
        """
        if len(triads) < 2:
            print(" *** make_udelm: At least two triads must be specified")
            return None

        base_ids = []
        for i in range(len(triads) - 1):
            base_id = self._fmlib.FmCreateUDE2(
                _convert_char(name),
                self._convert_id(triads[i], FmType.TRIAD),
                self._convert_id(triads[i + 1], FmType.TRIAD),
            )
            if base_id < 1:
                return None
            base_ids.append(base_id)

            if "alpha1" in kwargs or "alpha2" in kwargs:
                self._fmlib.FmStructDamp(
                    base_id,
                    _convert_real(kwargs.get("alpha1", 0.0)),
                    _convert_real(kwargs.get("alpha2", 0.0)),
                )

        if "tag" in kwargs:
            self.fm_tag_object(base_ids, _convert_char(kwargs["tag"]))

        return base_ids

    def make_assembly(self, name, objs=None):
        """
        Creates a sub-assembly and moves the specified objects into it.

        Parameters
        ----------
        name : str
            Description of the new sub-assembly
        objs : list of int or str or list of str, default=None
            List of base Ids or tags of the objects to put into the new sub-assembly

        Returns
        -------
        int
            Base Id of the created sub-assembly
        """

        obj_ids = []
        if isinstance(objs, (list, tuple)):
            for obj in objs:
                id_ = self._convert_id(obj, FmType.ALL, True)
                if isinstance(id_, list):
                    obj_ids.extend(id_)
                else:
                    obj_ids.append(id_)
        elif isinstance(objs, str):
            obj_ids = self._convert_id(objs, FmType.ALL, True)
        elif isinstance(objs, (int, c_int)):
            obj_ids = [objs]

        num_obj_, obj_ids_ = _convert_int_array(obj_ids)
        return self._fmlib.FmCreateAssembly(_convert_char(name), num_obj_, obj_ids_)

    def _move_object(self, base_id, kwargs):
        """
        Moves an existing object identified by its base Id.
        """
        dofs_ = _extract_dofs(kwargs)
        if dofs_ is None:
            return True  # No movement specified, silently ignore

        return self._fmlib.FmMoveObject(
            _convert_int(base_id),
            dofs_,
            _convert_int(kwargs.get("tra_ref", 0)),
            _convert_int(kwargs.get("rot_ref", 0)),
        )

    def _constrain_object(self, base_id, kwargs):
        """
        Constrains DOFs in an existing object identified by its base Id.
        """
        if kwargs is None:
            return True  # No constraints specified, silently ignore

        dof_tags = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "All"]
        for key, value in kwargs.items():
            if key in dof_tags:
                dof_ = c_int(dof_tags.index(key))
                con_ = _convert_int(value, FmDofStat)
                print(f"Constraining {key} in object [{base_id.value}]", value)
                if not self._fmlib.FmConstrainObject(base_id, dof_, con_):
                    return False

        return True

    def _dof_property(self, base_id, prt, kwargs):
        """
        Assigns a DOF property in an existing object identified by its base Id.
        """
        if kwargs is None:
            return True  # No properties specified, silently ignore

        pr_type_ = c_int(prt)
        dof_tags = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]
        for key, value in kwargs.items():
            if key in dof_tags:
                dof_ = c_int(dof_tags.index(key))
                if prt in (1, 4):
                    fcn_ = _convert_int(value, int_or_none=True)
                else:
                    fcn_ = None
                if fcn_ is None:  # A constant property value is specified
                    val_ = _convert_real(value)
                    fcn_ = c_int(0)
                else:  # The user Id of a general function is specified
                    val_ = c_double(0)
                print(f"{key} property {prt} in object [{base_id.value}]", value, fcn_)
                if not self._fmlib.FmDofProperty(base_id, dof_, pr_type_, val_, fcn_):
                    return False

        return True

    def edit_triad(self, obj_id, **kwargs):
        """
        Modifies an existing triad.

        Parameters
        ----------
        obj_id : int or str or list of int or list of str
            Base Id or tag of the triad(s) to modify

        kwargs : dict
            Keyword arguments containing the new triad attributes to assign.
            Currently, the following keywords are recognized:

            * Tx, Ty, Tz : Translation offset in global X, Y, and Z-direction
            * Rx, Ry, Rz : Euler angles (in degrees) for rotating the Triad
            * tra_ref : Base Id of object used as translation reference
            * rot_ref : Base Id of object used as rotation reference
            * mass : Additional mass (and inertia) on the Triad
            * mass_func : User Id or base Id (if negative) of mass scaling function

            * constraints : Dictionary with keywords for constraing DOFs, i.e.,
              ``constraints={"Tx" : <value>, ..., "Rz" : <value>}``, where
              ``<value>`` can be any of the enum values :class:`enums.FmDofStat`.
              Only the DOFs that should be changed need to be specified.
              You can also use "All" as key which implies all DOFs in the triad.

            * init_vel : Dictionary with keywords specifying initial velocities,
              i.e., ``init_vel={"Tx" : <u01>, ..., "Rz" : <u06>}``, where
              ``<u0i>`` is the initial velocity to be assigned local dof `i`.

            * load : Dictionary with keywords specifying constant DOF loads,
              i.e., ``load={"Tx" : <f1>, ..., "Rz" : <f6>}``, where
              ``<fi>`` is the user Id of a general function defining the load
              in local dof `i` if the type is int, otherwise it is taken as the
              constant load to be assigned local dof `i`.

            * motion : Dictionary with keywords specifying prescribed motions,
              i.e., ``motion={"Tx" : <u1>, ..., "Rz" : <u6>}``, where
              ``<ui>`` is the user Id of a general function defining the motion
              in local dof `i` if the type is int, otherwise it is taken as the
              constant motion to be prescribed in local dof `i`.

        Returns
        -------
        bool
            True if the input is valid, otherwise False
        """

        # Convert tag to base Id or list of base Id
        base_id = self._convert_ids(obj_id, FmType.TRIAD)

        retval = True
        if isinstance(base_id, list):  # Process a list of triads
            for bid in base_id:
                retval &= self.edit_triad(bid, **kwargs)
            return retval

        retval = self._move_object(base_id, kwargs)
        if not retval:
            print(" *** Failed to move triad", base_id.value)

        mass = kwargs.get("mass", None)
        if mass is not None:
            if isinstance(mass, tuple):
                n_, m_ = _convert_float_array(mass)
            else:
                n_ = c_int(1)
                m_ = (c_double * 1)(_convert_real(mass))
            mfunc_ = _convert_int(kwargs.get("mass_function", 0))
            if not self._fmlib.FmAddMass(base_id, n_, m_, mfunc_):
                print(" *** Failed to add mass on triad", base_id.value)
                retval = False

        if not self._constrain_object(base_id, kwargs.get("constraints", None)):
            print(" *** Failed to constrain triad", base_id.value)
            retval = False

        if not self._dof_property(base_id, 0, kwargs.get("init_vel", None)):
            print(" *** Failed to assign initial velocity to triad", base_id.value)
            retval = False

        load = kwargs.get("load", None)
        if load is None:
            load = kwargs.get("motion", None)
        if not self._dof_property(base_id, 1, load):
            print(" *** Failed to assign load/motion to triad", base_id.value)
            retval = False

        return retval

    def edit_part(self, obj_id, **kwargs):
        """
        Modifies an existing FE part.

        Parameters
        ----------
        obj_id : int or str or list of int or list of str
            Base Id or tag of the part(s) to modify

        kwargs : dict
            Keyword arguments containing the new part attributes to assign.
            Currently, the following keywords are recognized:

            * Tx, Ty, Tz : Translation offset in global X, Y, and Z-direction
            * Rx, Ry, Rz : Euler angles (in degrees) for rotating the Part
            * tra_ref : Base Id of object used as translation reference
            * rot_ref : Base Id of object used as rotation reference
            * alpha1 : Mass-proportional damping coefficient
            * alpha2 : Stiffness-proportional damping coefficient
            * component_modes : Number of component modes
            * consistent_mass : If True, use consistent mass (lumped is default)
            * recovery : Flag for activating FE part recovery during solve
              (0=off, 1=stress recovery, 2=strain gage recovery, 3=both)

        Returns
        -------
        bool
            True if the input is valid, otherwise False
        """

        # Convert tag to base Id or list of base Id
        base_id = self._convert_ids(obj_id, FmType.FEPART)

        retval = True
        if isinstance(base_id, list):  # Process a list of parts
            for bid in base_id:
                retval &= self.edit_part(bid, **kwargs)
            return retval

        retval = self._move_object(base_id, kwargs)

        if "alpha1" in kwargs or "alpha2" in kwargs:
            retval = retval and self._fmlib.FmStructDamp(
                base_id,
                _convert_real(kwargs.get("alpha1", 0.0)),
                _convert_real(kwargs.get("alpha2", 0.0)),
            )

        if "component_modes" in kwargs or "consistent_mass" in kwargs:
            retval = retval and self._fmlib.FmReduceOpts(
                base_id,
                _convert_int(kwargs.get("component_modes", 0)),
                _convert_bool(kwargs.get("consistent_mass", False)),
            )

        if "recovery" in kwargs:
            retval = retval and self._fmlib.FmRecoverOpts(
                base_id, _convert_int(kwargs["recovery"]), c_bool(False)
            )

        if not retval:
            print(" *** Failed to edit properties on part", base_id.value)

        return retval

    def edit_joint(self, obj_id, **kwargs):
        """
        Modifies an existing joint.

        Parameters
        ----------
        obj_id : int or str or list of int or list of str
            Base Id or tag of the joint(s) to modify

        kwargs : dict
            Keyword arguments containing the new joint attributes to assign.
            Currently, the following keywords are recognized:

            * Tx, Ty, Tz : Translation offset in global X, Y, and Z-direction
            * Rx, Ry, Rz : Euler angles (in degrees) for rotating the joint
            * tra_ref : Base Id of object used as translation reference
            * rot_ref : Base Id of object used as rotation reference

            * constraints : Dictionary with keywords for constraing DOFs, i.e.,
              ``constraints={"Tx" : <value>, ..., "Rz" : <value>}``, where
              ``<value>`` can be any of the enum values :class:`enums.FmDofStat`.
              Only the DOFs that should be changed need to be specified.
              You can also use "All" as key which implies all DOFs in the joint.

            * init_vel : Dictionary with keywords specifying initial velocities,
              i.e., ``init_vel={"Tx" : <u01>, ..., "Rz" : <u06>}``, where
              ``<u0i>`` is the initial velocity to be assigned local dof `i`.

            * load : Dictionary with keywords specifying constant DOF loads,
              i.e., ``load={"Tx" : <f1>, ..., "Rz" : <f6>}``, where
              ``<fi>`` is the user Id of a general function defining the load
              in local dof `i` if the type is int, otherwise it is taken as the
              constant load to be assigned local dof `i`.

            * motion : Dictionary with keywords specifying prescribed motions,
              i.e., ``motion={"Tx" : <u1>, ..., "Rz" : <u6>}``, where
              ``<ui>`` is the user Id of a general function defining the motion
              in local dof `i` if the type is int, otherwise it is taken as the
              constant motion to be prescribed in local dof `i`.

            * spring : Dictionary with keywords specifying spring stiffnesses,
              i.e., ``spring={"Tx" : <k1>, ..., "Rz" : <k6>}``, where
              ``<ki>`` is the constant stiffness to be assigned local dof `i`.

            * damper : Dictionary with keywords specifying damping coefficients,
              i.e., ``load={"Tx" : <c1>, ..., "Rz" : <c6>}``, where
              ``<ci>`` is the constant damping to be assigned local dof `i`.

            * length : Dictionary with keywords specifying stress-free lengths,
              i.e., ``length={"Tx" : <l1>, ..., "Rz" : <l6>}``, where
              ``<li>`` is the user Id of a general function defining the
              stress-free length in local dof `i`.

        Returns
        -------
        bool
            True if the input is valid, otherwise False
        """

        # Convert tag to base Id or list of base Id
        base_id = self._convert_ids(obj_id, FmType.JOINT)

        retval = True
        if isinstance(base_id, list):  # Process a list of joints
            for bid in base_id:
                retval &= self.edit_joint(bid, **kwargs)
            return retval

        retval = self._move_object(base_id, kwargs)
        if not retval:
            print(" *** Failed to move joint", base_id.value)

        if not self._constrain_object(base_id, kwargs.get("constraints", None)):
            print(" *** Failed to constrain joint", base_id.value)
            retval = False

        if not self._dof_property(base_id, 0, kwargs.get("init_vel", None)):
            print(" *** Failed to assign initial velocity to joint", base_id.value)
            retval = False

        load = kwargs.get("load", None)
        if load is None:
            load = kwargs.get("motion", None)
        if not self._dof_property(base_id, 1, load):
            print(" *** Failed to assign load/motion to joint", base_id.value)
            retval = False

        if not self._dof_property(base_id, 2, kwargs.get("spring", None)):
            print(" *** Failed to assign spring stiffness to joint", base_id.value)
            retval = False

        if not self._dof_property(base_id, 3, kwargs.get("damper", None)):
            print(" *** Failed to assign damping coefficient to joint", base_id.value)
            retval = False

        if not self._dof_property(base_id, 4, kwargs.get("length", None)):
            print(" *** Failed to assign stress-free lengths to joint", base_id.value)
            retval = False

        return retval

    def edit_function(self, func_id, obj, var=None, dof=None):
        """
        Modifies an existing function by changing its argument.

        Parameters
        ----------
        func_id : int
            User Id (or base Id, if negative) of the function to modify
        obj : int or str or (int, int) or (str, str)
            Base Id or tag of object(s) producing the response quantity to use
        var : FmVar, default=None=0
            Response quantity to use as function argument
        dof : FmDof, default=None=0
            Which component of `var` to be used if a multi-DOF quantity

        Returns
        -------
        bool
            True if the input is valid, otherwise False
        """

        if isinstance(obj, tuple):
            obj1_ = self._convert_id(obj[0])
            obj2_ = self._convert_id(obj[1])
        else:
            obj1_ = self._convert_id(obj)
            obj2_ = c_int(0)

        retval = self._fmlib.FmSetFunctionArg(
            _convert_int(func_id),
            _convert_int(var),
            _convert_int(dof),
            obj1_,
            obj2_,
        )

        if not retval:
            print(" *** Failed to change argument of Function", func_id)
            print("     to", obj, var, dof)
        return retval

    def solver_setup(self, **kwargs):
        """
        Setting up solver parameters, t_quasi, t_inc and t_end.
        Switching off initial equilibrium by setting t_quasi negative.

        Parameters
        ----------
        kwargs : dict
            Keyword arguments containing the solver settings to assign.
            Currently, the following keywords are recognized:

            * t_start : Start time
            * t_end : Stop time
            * t_inc : Time step size
            * t_quasi : Stop time for quasi-static simulation.
                        If equal to t_start, perform initial equilibrium analysis
                        before starting the dynamics time integration.
            * n_modes : If non-zero, perform eigenvalue analysis during the simulation,
                        and calculate this number of modes each time
            * e_inc: Time between each eigenvalue analysis
            * add_opt: Additional solver options
            * tol_ene: Convergence tolerance in energy norm
            * tol_dis: Convergence tolerance in displacement norm
            * tol_vel: Convergence tolerance in velocity norm
            * tol_res: Convergence tolerance in force residual norm
        """

        self.fm_solver_setup(
            kwargs.get("t_start", 0.0),
            kwargs.get("t_end", 1.0),
            kwargs.get("t_inc", 0.01),
            kwargs.get("t_quasi", 0.0),
            kwargs.get("n_modes", 0),
            kwargs.get("e_inc", 0.0),
            kwargs.get("add_opt", None),
        )

        if len(set(kwargs.keys()) & {"tol_ene", "tol_dis", "tol_vel", "tol_res"}) > 0:
            self.fm_solver_tol(
                kwargs.get("tol_ene", -1.0),
                kwargs.get("tol_dis", -1.0),
                kwargs.get("tol_vel", -1.0),
                kwargs.get("tol_res", -1.0),
            )
