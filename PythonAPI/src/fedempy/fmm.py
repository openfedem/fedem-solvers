# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Python wrapper for the Fedem Mechanism Model database library.
Used for convenience in order to hide native type convertions.
"""

from ctypes import byref, c_bool, c_char_p, c_double, c_int, cdll, create_string_buffer

from fedempy.enums import Enum, FmType


def _convert_char(arg):
    """
    Converts a string argument from Python to C type.
    """
    if arg is None:
        return c_char_p()

    if isinstance(arg, str):
        return c_char_p(arg.encode("utf-8"))

    if isinstance(arg, c_char_p):
        return arg

    raise TypeError(f"Expected {str} or {c_char_p}, got {type(arg)}.")


def _convert_int(arg, enum_type=None, int_or_none=False):
    """
    Converts an integer from Python to C type.
    If `enum_type` is specified, then `arg` needs to be an enum of that type.
    If `int_or_none` is True, then `None` is returned if `arg` is not an int,
    instead to raising an exception.
    """
    if arg is None:
        return c_int(0)

    if enum_type is not None:
        if isinstance(arg, enum_type):
            return c_int(arg.value)
        raise TypeError(f"Expected {enum_type}, got {type(arg)}.")

    if isinstance(arg, Enum):
        return c_int(arg.value)

    if isinstance(arg, int):
        return c_int(arg)

    if isinstance(arg, c_int):
        return arg

    if int_or_none:
        return None

    raise TypeError(f"Expected {Enum}, {int} or {c_int}, got {type(arg)}.")


def _convert_real(arg):
    """
    Converts a float from Python to C type.
    """
    if arg is None:
        return c_double(0)

    if isinstance(arg, (float, int)):
        return c_double(arg)

    if isinstance(arg, c_double):
        return arg

    raise TypeError(f"Expected {int}, {float} or {c_double}, got {type(arg)}.")


class FedemModel:
    """
    This class wraps the Fedem model file handling functionality.

    Parameters
    ----------
    lib_path : str
        Absolute path to the Fedem model database shared object library
    plugin1 : str, default=None
        Absolute path to user-defined element plugin (optional)
    plugin2 : str, default=None
        Absolute path to user-defined function plugin (optional)

    Methods
    -------
    fm_open:
        Opens the specified Fedem model file
    fm_save:
        Saves current model with updated results to file
    fm_new:
        Initializes an empty model
    fm_close:
        Closes current model
    fm_count:
        Counts mechanism objects of specified type
    fm_get_objects:
        Returns a list of base Ids for all objects of specified type
    fm_tag_object:
        Tags the specified object(s) in the model
    fm_get_pos:
        Returns the global position of the specified object
    fm_get_node:
        Returns the FE node number matching a spatial point
    fm_sync_part:
        Syncronizes an FE part with content on disk
    fm_write_reducer:
        Writes reducer input files for an FE part
    fm_write_solver:
        Writes solver input files
    fm_solver_setup:
        Defines some basis solver setup parameters
    fm_solver_tol:
        Defines the solver convergence tolerances
    def fm_get_func_tag:
        Returns the tag of an indexed external function
    """

    def __init__(self, lib_path, plugin1=None, plugin2=None):
        """
        Constructor.
        Initializes the internal datastructure of the shared object library.
        """
        self._fmlib = cdll.LoadLibrary(lib_path)

        self._fmlib.FmInit(_convert_char(plugin1), _convert_char(plugin2))
        self._fmlib.FmOpen.restype = c_bool
        self._fmlib.FmOpen.argtypes = [c_char_p]
        self._fmlib.FmSave.restype = c_bool
        self._fmlib.FmSave.argtypes = [c_char_p]
        self._fmlib.FmCount.restype = c_int
        self._fmlib.FmCount.argstype = [c_int]
        self._fmlib.FmGetObjects.restype = c_int
        self._fmlib.FmTagObjects.restype = c_int
        self._fmlib.FmGetPosition.restype = c_bool
        self._fmlib.FmGetNode.restype = c_int
        self._fmlib.FmGetFEModel.restype = c_bool
        self._fmlib.FmSync.restype = c_bool
        self._fmlib.FmSync.argstype = [c_int]
        self._fmlib.FmReduce.restype = c_bool
        self._fmlib.FmReduce.argtypes = [c_char_p, c_int]
        self._fmlib.FmSolve.restype = c_bool
        self._fmlib.FmSolve.argtypes = [c_char_p, c_bool, c_char_p, c_char_p]
        self._fmlib.FmGetFuncTag.restype = c_bool
        self._fmlib.FmGetFuncTag.argtypes = [c_int, c_char_p]

    def _convert_id(self, obj_id, obj_type=FmType.ALL, multiple=False):
        """
        Converts `obj_id`, which can be either an integer base Id or a string tag,
        into C type base Id. If `multiple` is True and `obj_id` is a tag,
        a list of base Ids is returned, even if only one (or none) item matches.
        """
        if not isinstance(obj_id, str):
            # The base Id is specified
            return _convert_int(obj_id)

        # The object tag is specified
        base_id = self.fm_get_objects(obj_type, obj_id)
        if multiple:
            return [c_int(bid) for bid in base_id]

        if len(base_id) == 1:
            return c_int(base_id[0])

        # Expected one object, but found none or many
        return c_int(0)

    def fm_open(self, fname):
        """
        This method opens the specified Fedem model file
        and reads its content into the internal datastructure.

        Parameters
        ----------
        fname : str
            Absolute path to the fmm-file to read

        Returns
        -------
        bool
            True on success, otherwise False
        """
        if fname is None:
            return False

        return self._fmlib.FmOpen(_convert_char(fname))

    def fm_save(self, fname=None):
        """
        This method saves the currently open Fedem model such that the model
        file is updated to reflect the current results data found on disk.

        Parameters
        ----------
        fname : str, default=None
            Absolute path to the fmm-file to write to.
            If not given, the current model file as specified by the last call
            to fm_open() or fm_new() is overwritten.

        Returns
        -------
        bool
            True on success, otherwise False
        """
        return self._fmlib.FmSave(_convert_char(fname))

    def fm_new(self, fname=None):
        """
        This method initializes an empty Fedem model.

        Parameters
        ----------
        fname : str, default=None
            Absolute path to the fmm-file to write to on next save.
            If not given, the name "untitled.fmm" will be used.

        Returns
        -------
        bool
            Always True
        """
        self._fmlib.FmNew(_convert_char(fname))

        return True

    def fm_close(self, remove_singletons=False):
        """
        This method closes the currently open Fedem model to clean up memory.

        Parameters
        ----------
        remove_singletons : bool, default=False
            If True, heap-allocated singelton objects are also released

        Returns
        -------
        bool
            Always True
        """

        if isinstance(remove_singletons, bool):
            self._fmlib.FmClose(c_bool(remove_singletons))
        elif isinstance(remove_singletons, c_bool):
            self._fmlib.FmClose(remove_singletons)
        else:
            self._fmlib.FmClose(c_bool(False))

        return True

    def fm_count(self, object_type=FmType.ALL):
        """
        This method returns the number of objects of the given type.

        Parameters
        ----------
        object_type : FmType, default=ALL
            Enum value identifying which object type to count the instances of

        Returns
        -------
        int
            Number of instances of the indicated object type
        """
        return self._fmlib.FmCount(_convert_int(object_type, FmType))

    def fm_get_objects(self, object_type=FmType.ALL, tag=None):
        """
        This method returns the base Id of all objects of the given type,
        and/or those with the specified tag.

        Parameters
        ----------
        object_type : FmType, default=ALL
            Enum value identifying which object type to return base Id of
        tag : str or list of str, default=None
            Return objects having this tag only, unless None

        Returns
        -------
        list
            Base Id list of all instances in the model of the specified type
        """
        base_id = []
        if isinstance(tag, list):  # Process a list of tags
            for atag in tag:
                base_id.extend(self.fm_get_objects(object_type, atag))
            return base_id

        num_objs = self.fm_count(object_type)
        base_id_ = (c_int * num_objs)()
        obj_typ_ = c_int(object_type.value)
        num_objs = self._fmlib.FmGetObjects(base_id_, obj_typ_, _convert_char(tag))
        base_id = [0] * num_objs
        base_id[:] = [base_id_[i] for i in range(num_objs)]

        return base_id

    def fm_tag_object(self, base_id, tag):
        """
        This method tags the specified object(s) identified by the base Id.

        Parameters
        ----------
        base_id : int or list of int
            Base Id of the object(s) to tag
        tag : str
            The tag to assign

        Returns
        -------
        int
            Number of tagged objects
        """
        if isinstance(base_id, list):
            num_obj_ = c_int(len(base_id))
            base_id_ = (c_int * len(base_id))()
            base_id_[:] = base_id
        else:
            num_obj_ = c_int(1)
            base_id_ = (c_int * 1)()
            base_id_[0] = base_id

        return self._fmlib.FmTagObjects(base_id_, num_obj_, _convert_char(tag))

    def fm_get_pos(self, object_id):
        """
        This method returns the global position of the specified object.

        Parameters
        ----------
        object_id : int or str
            Base Id or tag of the object to return the position for

        Returns
        -------
        (float, float, float)
            Global X-, Y-, and Z-coordinate
        """
        pos_ = (c_double * 3)()
        if not self._fmlib.FmGetPosition(self._convert_id(object_id), pos_):
            return None

        return (pos_[0], pos_[1], pos_[2])

    def fm_get_node(self, part_id, pos):
        """
        This method returns the FE node number matching a spatial point.

        Parameters
        ----------
        part_id : int
            Base Id or tag of the FE part that the node belongs to
        pos : (float, float, float)
            X-, Y-, and Z-coordinate of the node

        Returns
        -------
        int
            Node Id
        """
        pos_ = (c_double * 3)()
        pos_[:] = pos[0:3]
        return self._fmlib.FmGetNode(self._convert_id(part_id), pos_)

    def fm_get_femodel(self, base_id):
        """
        This method returns the FE data file for the specified Part.

        Parameters
        ----------
        base_id : int
            Base Id of the FE part to get the FE model for

        Returns
        -------
        str
            Absolute path to the FE data file
        int
            If less than -1: No file is associated with this part.
            If -1: This is a generic part with a visualization file.
            If 0: This is a FE part but not subjected to stress recovery.
            If 1: This FE part is subjected to stress recovery.
        """
        femfile = create_string_buffer(512)
        recover = c_int(0)
        if self._fmlib.FmGetFEModel(_convert_int(base_id), femfile, byref(recover)):
            return femfile.value.decode("utf-8"), recover.value

        return "", -2

    def fm_sync_part(self, base_id):
        """
        This method syncronizes the currently open Fedem model with contents
        on disk for the specified FE part. It is typically used after the
        reduction process for an FE part is finished.

        Parameters
        ----------
        base_id : int
            Base Id of the FE part to syncronize the model and disk content for

        Returns
        -------
        bool
            True on success, otherwise (e.g., if base_id is invalid) False
        """
        return self._fmlib.FmSync(_convert_int(base_id))

    def fm_write_reducer(self, base_id):
        """
        This method writes reducer input files for the specified FE part,
        unless that part is already reduced.

        Parameters
        ----------
        base_id : int
            Base Id of the FE part to generate reducer input files for

        Returns
        -------
        str
            Absolute path to the working directory of the reduction process
        """
        rdbdir = create_string_buffer(512)
        if self._fmlib.FmReduce(rdbdir, _convert_int(base_id)):
            return rdbdir.value.decode("utf-8")

        return ""

    def fm_write_solver(self, keep_old_res=False, ude=None, udf=None):
        """
        This method writes solver input files for the currently loaded model.

        Parameters
        ----------
        keep_old_res : bool, default=False
            Option to not overwrite any existing res-file in the RDB directory
        ude : str, default=None
            Absolute path to user-defined element plugin (optional)
        udf : str, default=None
            Absolute path to user-defined function plugin (optional)

        Returns
        -------
        str
            Absolute path to the working directory of the solver process
        """
        rdbdir = create_string_buffer(512)
        if self._fmlib.FmSolve(
            rdbdir, c_bool(keep_old_res), _convert_char(ude), _convert_char(udf)
        ):
            return rdbdir.value.decode("utf-8")

        return ""

    def fm_solver_setup(
        self,
        t_start=0.0,
        t_end=1.0,
        t_inc=0.01,
        t_quasi=0.0,
        n_modes=0,
        e_inc=0.0,
        add_opt=None,
    ):
        """
        This method (re)defines some basis solver setup parameters.

        Parameters
        ----------
        t_start : float, default=0
            Start time
        t_end : float, default=1
            Stop time
        t_inc : float, default=0.01
            Time step size
        t_quasi : float, default=0
            Stop time for quasi-static simulation.
            If equal to t_start, also perform initial equilibrium analysis
            before starting the time stepping.
        n_modes : int, default=0
            If non-zero, perform eigenvalue analysis during the simulation,
            and calculate this number of modes each time
        e_inc : float, default=0
            Time between each eigenvalue analysis
        add_opt : str, default=None
            Additional solver options
        """

        self._fmlib.FmSolveSetup(
            _convert_real(t_start),
            _convert_real(t_inc),
            _convert_real(t_end),
            _convert_real(t_quasi),
            _convert_real(e_inc),
            _convert_int(n_modes),
            _convert_char(add_opt),
        )

    def fm_solver_tol(self, tol_ene=-1.0, tol_dis=-1.0, tol_vel=-1.0, tol_res=-1.0):
        """
        This method (re)defines the solver convergence tolerances.

        Parameters
        ----------
        tol_ene : float, default=-1.0
            Energy norm tolerance
        tol_dis : float, default=-1.0
            Displacement norm tolerance
        tol_vel : float, default=-1.0
            Velocity norm tolerance
        tol_res : float, default=-1.0
            Residual force norm tolerance
        """

        self._fmlib.FmSolverTol(
            _convert_real(tol_ene),
            _convert_real(tol_dis),
            _convert_real(tol_vel),
            _convert_real(tol_res),
        )

    def fm_get_func_tag(self, channel):
        """
        This method returns the tag of an indexed external function.

        Parameters
        ----------
        channel : int
            Channel index of the external function

        Returns
        -------
        str
            The tag of the function, None if channel is out of range
        """

        tag = create_string_buffer(512)
        if self._fmlib.FmGetFuncTag(_convert_int(channel), tag):
            return tag.value.decode("utf-8")

        return None
