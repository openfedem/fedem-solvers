# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Python wrapper for the native fedem dynamics solver.
Used for convenience in order to hide native type convertions.
"""

from ctypes import byref, c_bool, c_char_p, c_double, c_int, cdll
from os import path

from numpy import empty, float64, int32, int64, ndarray
from progress.bar import Bar


class FedemProgressBar(Bar):
    """
    Progress bar for the dynamics solver.

    Parameters
    ----------
    solver : FedemSolver
        The dynamics solver object to report progress for
    """

    suffix = "%(percent).1f%% - %(elapsed)ds"

    def __init__(self, solver):
        """
        Constructor.
        """
        super().__init__("     Solving...")
        self._solver = solver
        self._tstart = solver.get_start_time()
        self._trange = solver.get_stop_time() - self._tstart

    @property
    def progress(self):
        """
        Override the progress property using the physical time of the simulation.
        """
        current_time = self._solver.get_current_time() - self._tstart
        return current_time / self._trange


class FedemException(Exception):
    """
    General exception-type for solver exceptions.

    Parameters
    ----------
    errmsg : str
        Error message to print.
    """

    def __init__(self, errmsg):
        """
        Constructor.
        """
        super().__init__({"Error": errmsg})


class FedemSolver:
    """
    This class mirrors the functionality of the fedem dynamics solver library
    (libfedem_solver_core.so on Linux, fedem_solver_core.dll on Windows).
    See also the header file ../../../src/vpmSolver/solverInterface.h

    In addition to the methods for conducting the simulation itself, the class
    also contains several methods for accessing and manipulating the linearized
    equation system, to facilitate use by external solution algorithms.

    Parameters
    ----------
    lib_path : str
        Absolute path to the solver shared object library
    solver_options : list of str, default=None
        List of command-line arguments passed to the solver
    use_internal_state : bool, default=False
        If True, internal state arrays are allocated (for microbatching)

    Methods
    -------
    solver_init:
        Processes the input and sets up necessary data structures
    restart_from_state:
        Re-initializes the mechanism objects with data from state array
    solve_window:
        Solves the problem for a time/load step window
    get_state_size:
        Returns the length of the state vector
    get_gauge_size:
        Returns the length of the initial strain gauge array
    get_transformation_state_size:
        Returns the length of the state vector holding transformation matrices
    get_part_deformation_state_size:
        Returns the length of the state vector holding deformation data
    get_part_stress_state_size:
        Returns the length of the state vector holding von Mises stress data
    save_state:
        Stores current solver state in the self.state_data array
    save_gauges:
        Stores initial gauge strains in the self.gauge_data array
    save_transformation_state:
        Stores current transformation state in provided core array
    save_part_state:
        Stores current deformation- and stress states in provided core arrays
    solve_next:
        Advances the solution one time/load step forward
    start_step:
        Starts a new time (or load) step
    solve_iteration:
        Solves current linearized equation system and updates state variables
    finish_step:
        Completes current time (or load) step
    solve_modes:
        Solves the eigenvalue problem at current time step
    solve_inverse:
        Solves the inverse problem at current time/load step
    solver_done:
        Closes down the model and cleans up heap memory and things on disk
    solver_close:
        Cleans up heap memory (singleton objects) on close
    run_all:
        Runs through the dynamics solver without any user intervention.
    set_ext_func:
        Assigns new value to an external function
    get_current_time:
        Returns the current physical time of the simulation
    get_next_time:
        Returns the physical time of the next step of the simulation
    get_start_time:
        Returns the start time of the simulation
    get_stop_time:
        Returns the stop time of the simulation
    get_function:
        Evaluates a general function in the model and returns its value
    get_functions:
        Evaluates several general functions in the model and returns their value
    get_function_ids:
        Returns a list of user Ids of tagged general functions
    get_equations:
        Returns the equation numbers associated with the DOFs of an object
    get_system_size:
        Returns the number of equations in the linearized system
    get_system_dofs:
        Returns the number of DOFs in the system
    get_newton_matrix:
        Returns current content of the system Newton matrix
    get_stiffness_matrix:
        Returns current content of the system stiffness matrix
    get_mass_matrix:
        Returns current content of the system mass matrix
    get_damping_matrix:
        Returns current content of the system damping matrix
    get_element_stiffness_matrix:
        Returns the content of a (beam) element stiffness matrix
    get_rhs_vector:
        Returns the content of the system right-hand-side vector
    get_external_force_vector:
        Returns the content of the external force vector
    set_rhs_vector:
        Replaces current content of the system right-hand-side vector
    add_rhs_vector:
        Updates the content of the system right-hand-side vector
    compute_strains_from_displ:
        Computes the strain tensor at gauges for given displacement field
    get_current_strains:
        Returns the current strain tensor for the specified gauges
    compute_rel_dist_from_displ:
        Computes relative distance at sensors for given displacement field
    compute_int_forces_from_displ:
        Computes beam section forces in triads for given displacement field
    compute_spring_var_from_displ:
        Computes one of the spring variables for given displacement field
    get_joint_spring_stiffness:
        Returns current joint spring stiffness coefficient(s)

    """

    def __init__(self, lib_path, solver_options=None, use_internal_state=False):
        """
        Constructor.
        Optionally initializes the solver itself if solver_options is given.
        """
        # load the solver library
        self._solver = cdll.LoadLibrary(lib_path)

        # set up return type for functions in the solver library
        self._solver.solverInit.restype = c_int
        self._solver.restartFromState.restype = c_int
        self._solver.solveWindow.restype = c_bool
        self._solver.haveResults.restype = c_int
        self._solver.getStateSize.restype = c_int
        self._solver.getTransformationStateSize.restype = c_int
        self._solver.getPartDeformationStateSize.restype = c_int
        self._solver.getPartStressStateSize.restype = c_int
        self._solver.getGagesSize.restype = c_int
        self._solver.saveState.restype = c_bool
        self._solver.saveTransformationState.restype = c_bool
        self._solver.savePartDeformationState.restype = c_bool
        self._solver.savePartStressState.restype = c_bool
        self._solver.saveGages.restype = c_bool
        self._solver.solveNext.restype = c_bool
        self._solver.startStep.restype = c_bool
        self._solver.solveIteration.restype = c_bool
        self._solver.solveEigenModes.restype = c_bool
        self._solver.solveInverse.restype = c_bool
        self._solver.solverDone.restype = c_int
        self._solver.setExtFunc.restype = c_bool
        self._solver.getTime.restype = c_double
        self._solver.evalFunc.restype = c_double
        self._solver.getEquations.restype = c_int
        self._solver.getStateVar.restype = c_int
        self._solver.getSystemSize.restype = c_int
        self._solver.getSystemMatrix.restype = c_bool
        self._solver.getElementStiffnessMatrix.restype = c_bool
        self._solver.getRhsVector.restype = c_bool
        self._solver.setRhsVector.restype = c_bool
        self._solver.addRhsVector.restype = c_bool
        self._solver.getBeamForcesFromDisp.restype = c_bool
        self._solver.getStrainsFromDisp.restype = c_bool
        self._solver.getRelDisp.restype = c_bool
        self._solver.getRespVars.restype = c_bool
        self._solver.getJointSprCoeff.restype = c_bool

        # initialize error flag
        self.ierr = c_int(-999)

        # initialize the internal state arrays
        self.state_size = c_int(0) if use_internal_state else c_int(-1)
        self.gauge_size = c_int(0)
        self.state_data = None
        self.gauge_data = None

        # initialize the fedem solver
        status = self.solver_init(solver_options)
        if status < 0:
            raise FedemException(
                f"Initialization failure ({status}). "
                + "Check the fedem_solver.res file for error messages."
            )

    def __check_error(self, func):
        """
        Checks that the internal error flag is zero, to prevent invoking another
        solver method if an error condition has occurred in a previous call.
        """
        if self.ierr.value != 0:
            raise FedemException(
                func + f"() cannot be called due to previous error ({self.ierr.value})."
            )

    @staticmethod
    def _convert_c_double(arg, default_value=None):
        """
        Converts a float from Python to C type.
        """
        if arg is None and default_value is not None:
            return c_double(default_value)
        if type(arg) in (float, int, float64, int32, int64):
            return c_double(arg)
        if isinstance(arg, c_double):
            return arg

        raise TypeError(
            f"Expected {float}, {int}, {float64}, {int32} or {int64}, got {type(arg)}."
        )

    @staticmethod
    def _convert_c_int(arg):
        """
        Converts an integer from Python to C type.
        """
        if type(arg) in (int, int32, int64):
            return c_int(arg)
        if isinstance(arg, c_int):
            return arg

        raise TypeError(f"Expected {int}, {int32} or {int64}, got {type(arg)}.")

    @staticmethod
    def _convert_c_char(arg):
        """
        Converts a character string from Python to C type.
        """

        if arg is None:
            return c_char_p(0)
        if isinstance(arg, str):
            return c_char_p(arg.encode("utf-8"))
        if isinstance(arg, c_char_p):
            return arg

        raise TypeError(f"Expected {str} or {c_char_p}, got {type(arg)}.")

    @staticmethod
    def _convert_c_double_array(arg, allow_none=False, ndiv=None):
        """
        Converts an array of floats from Python to C type.
        """
        if arg is None and allow_none:
            return c_int(0), None

        if type(arg) in (list, ndarray):
            argc = len(arg)
            argv = (c_double * argc)()
            argv[:] = arg
            if ndiv is None:
                return c_int(argc), argv
            return c_int(argc // ndiv), argv

        if hasattr(arg, "_length_") and getattr(arg, "_type_", None) is c_double:
            return c_int(len(arg)), arg

        raise TypeError(f"Expected {list}, {ndarray} or {c_double}, got {type(arg)}.")

    @staticmethod
    def _convert_c_int_array(arg, allow_none=False):
        """
        Converts an array of integers from Python to C type.
        """
        if arg is None and allow_none:
            return c_int(0), None

        if type(arg) in (list, ndarray):
            argc = len(arg)
            argv = (c_int * argc)()
            argv[:] = arg
            return c_int(argc), argv

        if hasattr(arg, "_length_") and getattr(arg, "_type_", None) is c_int:
            return c_int(len(arg)), arg

        raise TypeError(f"Expected {list}, {ndarray} or {c_int}, got {type(arg)}.")

    @staticmethod
    def _convert_c_char_array(arg):
        """
        Converts an array of strings from Python to C type.
        """
        if arg is None:
            return c_int(0), None

        if isinstance(arg, list):
            argc = len(arg)
            argv = (c_char_p * argc)()
            argv[:] = [aval.encode("utf-8") for aval in arg]
            return c_int(argc), argv

        if hasattr(arg, "_length_") and getattr(arg, "_type_", None) is c_char_p:
            return c_int(len(arg)), arg

        raise TypeError(f"Expected {list} or {c_char_p}, got {type(arg)}.")

    def solver_init(
        self, options, fsi=None, state_data=None, gauge_data=None, extf_input=None
    ):
        """
        This method processes the input and sets up necessary data structures
        prior to the time integration loop. It also performs the license checks.

        See the Fedem R8.0 Users Guide, Appendix C.3 for a complete list of all
        command-line arguments that may be specified and their default values.

        If the argument fsi is None, the model is assumed defined in the file
        specified via the command-line argument -fsifile instead.

        Parameters
        ----------
        options : list of str
            List of command-line arguments passed to the solver
        fsi : str, default=None
            Content of the solver input file describing the model
        state_data : list of float, default=None
            Complete state vector to restart simulation from
        gauge_data : list of float, default=None
            Initial strain gauge values for restart
        extf_input : list of float, default=None
            Initial external function values, for initial equilibrium iterations

        Returns
        -------
        int
            A negative value indicates an error, otherwise success.
            A positive value indicates that initial equilibrium iterations
            was performed with external function values from `extf_input`
        """
        if options is None:
            return 0  # do nothing if no solver options specified

        # The first option is expected to contain the path of the executable,
        # so insert a dummy name here if the first option starts with '-'
        if isinstance(options, list) and options[0][0] == "-":
            options.insert(0, "fedem_solver")  # arbitrary name (not used)
        argc_, argv_ = self._convert_c_char_array(options)

        cfsi_ = self._convert_c_char(fsi)

        ndat_, sdat_ = self._convert_c_double_array(state_data, True)
        ngda_, gdat_ = self._convert_c_double_array(gauge_data, True)
        nxin_, xinp_ = self._convert_c_double_array(extf_input, True)

        status = self._solver.solverInit(
            argc_, argv_, cfsi_, sdat_, ndat_, gdat_, ngda_, xinp_, byref(nxin_)
        )
        if status < 0:
            self.ierr = c_int(status)
            return status  # initialization failure

        self.ierr = c_int(0)
        if self.state_size.value < 0:
            return nxin_.value  # not using internal state arrays (no restart)

        # Initialize state array sizes
        self.state_size = c_int(self.get_state_size())
        self.gauge_size = c_int(self.get_gauge_size())

        # Initialize the state arrays
        if self.state_data is None:
            self.state_data = (c_double * self.state_size.value)()
        elif len(self.state_data) != self.state_size.value:
            raise FedemException(f"Invalid valid array size {len(self.state_data)}")
        if self.gauge_size.value > 0:
            if self.gauge_data is None:
                self.gauge_data = (c_double * self.gauge_size.value)()
            elif len(self.gauge_data) != self.gauge_size.value:
                raise FedemException(f"Invalid valid array size {len(self.gauge_data)}")

        return nxin_.value

    def restart_from_state(self, state_data, write_to_rdb=2):
        """
        This method re-initializes the mechanism objects with data from the
        provided state array, such that the simulation can continue from there.

        Parameters
        ----------
        state_data : list of float
            Complete state vector to restart simulation from
        write_to_rdb : int, default=2
            | Flag for saving response variables to results database,
            | = 0 : No results saving,
            | = 1 : Append results to already opened results database,
            | = 2 : Increment the results database and write to new files

        Returns
        -------
        int
            Zero on success, a negative value indicates some error
        """
        self.__check_error("restart_from_state")

        ndat_, sdat_ = self._convert_c_double_array(state_data)

        if write_to_rdb in (0, 1, 2):
            write_to_rdb_ = c_int(write_to_rdb)
        else:
            raise TypeError(f"Expected integer value in (0, 1, 2), got {write_to_rdb}.")

        return self._solver.restartFromState(sdat_, ndat_, write_to_rdb_)

    def solve_window(self, n_step, inputs=None, f_out=None, xtimes=None):
        """
        This method solves the problem for a time/load step window,
        with given values for the external functions, and extraction
        of results from another set of general functions in the model.

        A non-zero value on self.ierr on exit indicates that an error condition
        that will require the simulation to terminate has occurred.

        Parameters
        ----------
        n_step : int
            Number of time/load steps to solve for from current state
        inputs : list of float, default=None
            List of input sensor values for each time step.
            The length of this list must be equal to `n_step` times
            the number of input sensors.
        f_out : list of int, default=None
            List of user Ids identifying the output sensors in the model
        xtimes : list of float, default=None
            List of times associated with the inputs.
            The length of this list must be equal to `n_step`.

        Returns
        -------
        list of float
            Output sensor values for each time step
        bool
            Always True, unless the end of the simulation has been reached
        """
        self.__check_error("solve_window")

        n_inc_, xtimes_ = self._convert_c_double_array(xtimes, True)
        n_inp_, inputs_ = self._convert_c_double_array(inputs, True, n_step)
        n_out_, f_out_ = self._convert_c_int_array(f_out, True)
        outputs_ = (c_double * (n_step * len(f_out)))() if f_out else None

        not_done = self._solver.solveWindow(
            self._convert_c_int(n_step),
            n_inc_,
            n_inp_,
            n_out_,
            f_out_,
            xtimes_,
            inputs_,
            outputs_,
            self.state_size,
            self.state_data,
            byref(self.ierr),
        )

        if self.ierr.value != 0:
            success = False
        elif self.gauge_size.value > 0:
            success = self._solver.saveGages(self.gauge_data, self.gauge_size)
        else:
            success = True

        if f_out is None:
            outputs = None
        else:
            outputs = [0.0] * (len(f_out) * n_step)
            outputs[:] = outputs_

        return outputs, not_done and success

    def have_results(self):
        """
        Utility returning whether current time step have results to be saved.
        """
        return self._solver.haveResults()

    def get_state_size(self):
        """
        Utility returning the required length of the state vector
        which is used when restarting a simulation from an in-core array.
        """
        return self._solver.getStateSize()

    def get_gauge_size(self):
        """
        Utility returning the required size of the initial strain gauge array
        which is used when restarting a simulation from an in-core array.
        Returns 0 if the model does not contain any strain gauges.
        """
        return self._solver.getGagesSize()

    def get_transformation_state_size(self):
        """
        Utility returning the required length of the vector
        which stores the transformation matrices (rotation and translation)
        for Triads, Parts and Beams.

        | The size/length of the state vector is:
        | 3 +
        | (number of Triads) * 14 +
        | (number of Parts) * 14
        | (number of Beams) * 14
        """
        return self._solver.getTransformationStateSize()

    def get_part_deformation_state_size(self, base_id):
        """
        Utility returning the required length of the state vector
        which stores deformation data for the FE Part with the given base Id.
        Returns -1 if the specified Part does not exist.

        | The size/length of the state vector is:
        | (number of nodal points in the FE Part) * 3
        """
        bid_ = self._convert_c_int(base_id)
        return self._solver.getPartDeformationStateSize(bid_)

    def get_part_stress_state_size(self, base_id):
        """
        Utility returning the required length of the state vector
        which stores von Mises stresses for the FE Part with the given base Id.
        Returns -1 if the specified Part does not exist, and 0 if the specified
        Part does not contain any (shell or solid) elements with stresses.
        """
        bid_ = self._convert_c_int(base_id)
        return self._solver.getPartStressStateSize(bid_)

    def save_state(self):
        """
        This method stores current solver state in the self.state_data array.

        Returns
        -------
        bool
            Always True, unless the self.state_data array is too small
        """
        return self._solver.saveState(self.state_data, self.state_size)

    def save_gauges(self):
        """
        This method stores initial gauge strains in the self.gauge_data array.

        Returns
        -------
        bool
            Always True, unless the self.gauge_data array is too small
        """
        return self._solver.saveGages(self.gauge_data, self.gauge_size)

    def save_transformation_state(self, state_data):
        """
        This method stores current transformation state for Triads, Parts
        and Beams in the provided core array.

        The transformation state data is on the format:

        | [step number]
        | [current time]
        | [current time increment]
        | for each non-fixed Triad:
        |   1
        |   [rotMatrix column 1]
        |   [rotMatrix column 2]
        |   [rotMatrix column 3]
        |   [translation vector]
        | for each Part and Beam:
        |   2
        |   [rotMatrix column 1]
        |   [rotMatrix column 2]
        |   [rotMatrix column 3]
        |   [translation vector]

        Parameters
        ----------
        state_data : list of c_double
            Array to fill with transformation data

        Returns
        -------
        bool
            Always True, unless the state_data array is too small
        """
        ndat = c_int(len(state_data))
        return self._solver.saveTransformationState(state_data, ndat)

    def save_part_state(self, base_id, def_state, str_state):
        """
        This method stores current deformation- and stress states
        for the specified FE Part in the provided core arrays.

        Parameters
        ----------
        base_id : int
            Base Id of the FE Part to save state for
        def_state : list of c_double
            Array to fill with deformation data
        str_state : list of c_couble
            Array to fill with stress data

        Returns
        -------
        bool
            Always True, unless one or both of the state arrays are too small
        """
        bid_ = self._convert_c_int(base_id)
        ndef = c_int(len(def_state))
        if not self._solver.savePartDeformationState(bid_, def_state, ndef):
            return False

        nstr = c_int(len(str_state))
        return self._solver.savePartStressState(bid_, str_state, nstr)

    def solve_next(self, inp=None, inp_def=None, out_def=None, time_next=None):
        """
        This method advances the solution one time/load step forward.
        The self.ierr variable has the value zero on a successful computation.
        A non-zero value indicates some error that requires the simulation
        to terminate.

        Parameters
        ----------
        inp : list of float, default=None
            Input function values
        inp_def : list of int, default=None
            External function Ids of the functions to assign values
        out_def : list of int, default=None
            User Ids of the functions to evaluate the response for
        time_next : float, default=None
            Time of next step, to override time step size defined in the model

        Returns
        -------
        list of float, only if out_def is specified
            Evaluated response variables
        bool
            Always True, unless current time/load step failed to converge,
            or the end time of the simulation has been reached
        """
        self.__check_error("solve_next")
        if time_next is None:
            success = True
        else:
            success = self._solver.setTime(self._convert_c_double(time_next))

        if inp is not None:
            for i, val in enumerate(inp):
                if val is not None:
                    if inp_def is None:
                        self.set_ext_func(i + 1, val)
                    elif i < len(inp_def):
                        self.set_ext_func(inp_def[i], val)

        success &= self._solver.solveNext(byref(self.ierr))
        if out_def is None:
            return success

        return self.get_functions(out_def), success

    def start_step(self, time_next=None):
        """
        This method starts a new time (or load) step, by calculating the
        predicted response, the coefficient matrix and right-hand-side vector
        of the first nonlinear iteration. It has to be followed up by
        a series of solve_iteration calls in order to continue the simulation,
        but the linear equation system can be manipulated in between.
        The self.ierr variable has the value zero on a successful computation.
        A non-zero value indicates some error that requires the simulation
        to terminate.

        Parameters
        ----------
        time_next : float, default=None
            Time of next step, to override time step size defined in the model

        Returns
        -------
        bool
            Always True, unless the simulation has to stop due to some error,
            or the end time of the simulation has been reached
        """
        self.__check_error("start_step")
        if time_next is None:
            success = True
        else:
            success = self._solver.setTime(self._convert_c_double(time_next))

        return self._solver.startStep(byref(self.ierr)) if success else False

    def solve_iteration(self):
        """
        This method solves the current linearized equation system and updates
        all state variables. Then it assembles the linearized system of
        equations for next iteration, unless convergence has been reached.
        The self.ierr variable has the value zero on a successful computation.
        A non-zero value indicates some error that requires the simulation
        to terminate.

        Returns
        -------
        bool
            Always True, unless the simulation has to stop due to some error,
            or the current time/load step has converged.
        """
        self.__check_error("solve_iteration")
        return self._solver.solveIteration(byref(self.ierr), c_bool(False))

    def finish_step(self):
        """
        This method completes current time (or load) step, by iterating the
        linearized equation system until convergence is achieved.
        The self.ierr variable has the value zero on a successful computation.
        A non-zero value indicates some error that requires the simulation
        to terminate.

        Returns
        -------
        bool
            Always True, unless current time/load step failed to converge,
            or the end time of the simulation has been reached
        """
        self.__check_error("finish_step")
        return self._solver.solveIteration(byref(self.ierr), c_bool(True))

    def solve_modes(self, n_modes, dof_order=False, use_lapack=0):
        """
        This method solves the eigenvalue problem at current time step,
        and returns the computed eigenvalues and associated eigenvectors.
        If an error condition that requires the simulation to terminate occurs,
        the self.ierr variable is assigned a negative value.

        Parameters
        ----------
        n_modes : int
            Number of eigenmodes to calculate
        dof_order : bool, default=False
            If True, the eigenvectors are returned in DOF-order
            instead of equation order which is the default
        use_lapack : int, default=0
            Flag usage of LAPACK eigensolvers (0=No, 1=DSYGVX, 2=DGGEVX)

        Returns
        -------
        list of float
            The computed eigenvalues
        list of list of float
            The computed eigenvectors
        bool
            Always True, unless the computation failed
        """
        self.__check_error("solve_modes")

        dim = self.get_system_dofs() if dof_order else self.get_system_size()
        n_mod_ = self._convert_c_int(n_modes)
        e_val_ = (c_double * n_modes)()
        e_vec_ = (c_double * (dim * n_modes))()
        doford = c_bool(dof_order)
        lapack = c_int(use_lapack)
        success = self._solver.solveEigenModes(
            n_mod_, e_val_, e_vec_, doford, lapack, byref(self.ierr)
        )
        if self.ierr.value < 0:
            return None, None, success
        if self.ierr.value > 0:
            self.ierr = c_int(0)
            return None, None, success

        e_val = [0.0] * n_modes
        e_val[:] = e_val_
        e_vec = [e_vec_[dim * i : dim * (i + 1)] for i in range(n_modes)]

        return e_val, e_vec, success

    def solve_inverse(self, x_val, x_def, g_def, out_def=None):
        """
        This method solves the inverse problem at current time/load step,
        assuming small deformations only (linear response).
        If an error condition that requires the simulation to terminate occurs,
        the self.ierr variable is assigned a non-zero value.

        Parameters
        ----------
        x_val : list of float
            Specified displacement values at a set of degrees of freedom
        x_def : list of int
            Equation numbers for the specified displacement values
        g_def : list of int
            Equation numbers for the DOFs with unknown external forces
        out_def : list of int, default=None
            User Ids of the functions to evaluate the response for

        Returns
        -------
        list of float, only if out_def is specified
            Evaluated response variables
        bool
            Always True, unless the simulation has to stop due to some error,
            or the end of the simulation has been reached
        """
        self.__check_error("solve_inverse")

        nval_, xval_ = self._convert_c_double_array(x_val)
        ndis_, xeqs_ = self._convert_c_int_array(x_def)
        nfrs_, feqs_ = self._convert_c_int_array(g_def)
        if nval_.value < ndis_.value:
            raise FedemException(f"Array x_val is too small ({nval_.value}).")

        success = self._solver.solveInverse(
            xval_, xeqs_, feqs_, ndis_, nfrs_, byref(self.ierr)
        )
        if out_def is None:
            return success

        return self.get_functions(out_def), success

    def solver_done(self, remove_singletons=None, print_res=False):
        """
        This method should be used when the time/load step loop is finished.
        It closes down the model and cleans up heap memory and things on disk.

        Parameters
        ----------
        remove_singletons : bool, default=None
            If True or None, heap-allocated singelton objects are also released
        print_res : bool, default=False
            If True, the res-file content is printed to console

        Returns
        -------
        int
            Zero on success, non-zero values indicates errors.
        """
        if remove_singletons is None:
            _rsflag = c_bool(True)
        elif isinstance(remove_singletons, bool):
            _rsflag = c_bool(remove_singletons)
        elif isinstance(remove_singletons, c_bool):
            _rsflag = remove_singletons
        else:
            _rsflag = c_bool(True)

        status = self._solver.solverDone(_rsflag)
        if print_res and path.isfile("./fedem_solver.res"):
            with open("./fedem_solver.res", "r") as resf:
                print_line = 1
                for count, line in enumerate(resf):
                    if line.find("Error :") == 0:
                        print_line = 2
                    elif count > 100 and print_line == 1:
                        print_line = 0
                        print("     . . .")
                    if print_line > 0:
                        print("%8d %s" % (count + 1, line.rstrip()))
                print("#### End of file fedem_solver.res", flush=True)
        return status

    def solver_close(self):
        """
        This method needs to be used if solver_done() was invoked
        with its `remove_singletons` argument set to False.
        It will delete those singleton objects here instead.
        """
        self._solver.solverClose()

    def run_all(self, options):
        """
        This method runs the dynamics solver with given command-line `options`,
        without any user intervention.
        """
        status = self.solver_init(options)
        if status < 0:
            print(f" *** Solver initialization failure ({status}).", flush=True)
            print("     Check the fedem_solver.res file for error messages.")
            return status

        with FedemProgressBar(self) as pbar:
            pbar.next()
            while self.solve_next():
                pbar.next()
            if self.ierr.value == 0:
                pbar.next()

        if self.solver_done() == 0 and self.ierr.value == 0:
            print("     Time step loop OK")
        else:
            print(" *** Dynamics solver failed", self.ierr.value, flush=True)

        return self.ierr.value

    def set_ext_func(self, func_id, value=None):
        """
        This method may be used prior to the solve_next call, to assign a
        sensor value from a physical twin to the specified actuator or load in
        the model, identified by the argument func_id (external function Id).

        Parameters
        ----------
        func_id : int
            Id of the external function to be assigned new value
        value : float, default=None
            The value to be assigned

        Returns
        -------
        bool
            Always True, unless the specified function is not found
        """
        func_id_ = self._convert_c_int(func_id)
        f_value_ = self._convert_c_double(value, 0)
        return self._solver.setExtFunc(func_id_, f_value_)

    def get_current_time(self):
        """
        Utility returning the current physical time of the simulation.
        The self.ierr variable is not touched.
        """
        return self._solver.getTime(c_int(0), byref(self.ierr))

    def get_next_time(self):
        """
        Utility returning the physical time of the next step of the simulation.
        If the time step size is defined by a general function that could not be
        evaluated, the self.ierr variable is decremented.
        Otherwise, it is not touched.
        """
        return self._solver.getTime(c_int(1), byref(self.ierr))

    def get_start_time(self):
        """
        Utility returning the start time of the simulation.
        The self.ierr variable is not touched.
        """
        return self._solver.getTime(c_int(2), byref(self.ierr))

    def get_stop_time(self):
        """
        Utility returning the stop time of the simulation.
        The self.ierr variable is not touched.
        """
        return self._solver.getTime(c_int(3), byref(self.ierr))

    def check_times(self, xtimes, use_times=True):
        """
        Utility checking that a given time series starts with the current solver time.
        """
        if not use_times or type(xtimes) not in (list, ndarray):
            return None

        mismatch = xtimes[0] - self.get_next_time()
        if abs(mismatch) > 0.001:
            print(f"  ** Warning: Mismatch ({mismatch}) between the first input time")
            print(f"     and expected time of the first step ({self.get_next_time()}).")
            print("     This may cause problems and/or incorrect simulation results.")

        return xtimes

    def get_function(self, uid=0, tag=None, arg=None):
        """
        Utility evaluating a general function in the model, identified by the
        specified user Id `uid` or `tag`, and with the function argument `arg`.
        If `arg` is `None` or contains a negative value, the function is instead
        evaluated for current state of the sensor object that is defined as
        the function argument. The same is done if a `tag` is specified.
        If the specified function could not be evaluated, the self.ierr variable
        is decremented. Otherwise, it is not touched.
        """
        uid_ = self._convert_c_int(uid)
        tag_ = self._convert_c_char(tag)
        arg_ = self._convert_c_double(arg, -1)
        return self._solver.evalFunc(uid_, tag_, arg_, byref(self.ierr))

    def get_functions(self, uids):
        """
        Utility evaluating a list of general functions for current state.
        The `uids` argument can be a list of either the user Ids of the
        functions to be evaluated, or their corresponding objects tags.
        If the specified functions could not be evaluated, the self.ierr variable
        is decremented for each problem encountered. Otherwise, it is not touched.
        """
        out = [0.0] * len(uids)
        for i, uid in enumerate(uids):
            if isinstance(uid, str):
                out[i] = self.get_function(tag=uid)
            else:
                out[i] = self.get_function(uid)

        return out

    def get_function_ids(self, tags):
        """
        Utility returning a list of user Ids of tagged general functions.
        """
        if tags is None:
            return None
        if isinstance(tags, str):
            return self._solver.getFuncId(self._convert_c_char(tags))
        if not isinstance(tags, list):
            return None

        uids = [0] * len(tags)
        for i, tag in enumerate(tags):
            if isinstance(tag, str):
                uids[i] = self._solver.getFuncId(self._convert_c_char(tag))
            else:
                uids[i] = self._convert_c_int(tag)

        return uids

    def get_equations(self, bid):
        """
        Utility returning the equation numbers for the DOFs of the object
        with the specified base Id (bid).
        """
        bid_ = self._convert_c_int(bid)
        eqns = (c_int * 6)()
        num_eqn = self._solver.getEquations(bid_, eqns)
        meqn = [0] * num_eqn
        meqn[:] = eqns

        return meqn

    def get_system_size(self):
        """
        Utility returning the dimension (number of equations) of the system.
        """
        return self._solver.getSystemSize(c_bool(False))

    def get_system_dofs(self):
        """
        Utiloty returning the total number of DOFs of the system.
        """
        return self._solver.getSystemSize(c_bool(True))

    def __get_system_matrix(self, i_mat):
        """
        Utility returning a system matrix.
        """
        dim = self.get_system_size()
        matrix = (c_double * (dim * dim))()
        success = self._solver.getSystemMatrix(matrix, c_int(i_mat))
        n_mat = empty([dim, dim])

        for i in range(dim):
            for j in range(dim):
                n_mat[i, j] = matrix[i + j * dim]

        return n_mat, success

    def get_newton_matrix(self):
        """
        Utility returning current content of the system Newton matrix.
        """
        return self.__get_system_matrix(0)

    def get_stiffness_matrix(self):
        """
        Utility returning current content of the system stiffness matrix.
        """
        return self.__get_system_matrix(1)

    def get_mass_matrix(self):
        """
        Utility returning current content of the system mass matrix.
        """
        return self.__get_system_matrix(2)

    def get_damping_matrix(self):
        """
        Utility returning current content of the system damping matrix.
        """
        return self.__get_system_matrix(3)

    def get_element_stiffness_matrix(self, bid):
        """
        Utility returning the initial content of an element stiffness matrix.
        Use this method for beam elements only, 6 dof per node (total 12).
        """
        dim = 12
        matrix = (c_double * (dim * dim))()
        bid_ = self._convert_c_int(bid)
        success = self._solver.getElementStiffnessMatrix(matrix, bid_)
        es_mat = empty([dim, dim])

        for i in range(dim):
            for j in range(dim):
                es_mat[i, j] = matrix[i + j * dim]

        return es_mat, success

    def __get_rhs_vector(self, i_vec):
        """
        Utility returning a system right-hand-side vector.
        """
        dim = self.get_system_size()
        vec = (c_double * dim)()
        success = self._solver.getRhsVector(vec, c_int(i_vec))
        r_vec = empty(dim)

        for i in range(dim):
            r_vec[i] = vec[i]

        return r_vec, success

    def get_rhs_vector(self):
        """
        Utility returning current content of the system right-hand-side vector.
        """
        return self.__get_rhs_vector(0)

    def get_external_force_vector(self):
        """
        Utility returning current content of the external force vector.
        """
        return self.__get_rhs_vector(5)

    def set_rhs_vector(self, r_vec):
        """
        Utility replacing current content of the system right-hand-side vector.
        """
        return self._solver.setRhsVector(self._convert_c_double_array(r_vec)[1])

    def add_rhs_vector(self, r_vec):
        """
        Utility updating the content of the system right-hand-side vector.
        """
        return self._solver.addRhsVector(self._convert_c_double_array(r_vec)[1])

    def compute_strains_from_displ(self, disp, gauge_ids):
        """
        This method computes the strain tensor at gauges for the given
        displacement field, or from current state if no displacements provided.

        Parameters
        ----------
        disp : list of float
            Displacement field to evaluate strains from
        gauge_ids : list of int
            Array with gauge identification numbers

        Returns
        -------
        list of float
            Strain components
        bool
            Always True, unless the calculation failed
        """

        n_dof_, disp_ = self._convert_c_double_array(disp, True)

        ngauge = len(gauge_ids) // 2
        ngauge_ = c_int(ngauge)
        gauge_ids_ = (c_int * ngauge)()
        for i in range(ngauge):
            gauge_ids_[i] = self._convert_c_int(gauge_ids[2 * i])

        eps_ = (c_double * (3 * ngauge))()
        success = self._solver.getStrainsFromDisp(
            disp_, gauge_ids_, eps_, n_dof_, ngauge_
        )

        eps = []
        for i in range(ngauge):
            tc_idx = gauge_ids[2 * i + 1]
            if tc_idx in (0, 1, 2):  # Tensorial component index
                eps.append(eps_[3 * i + tc_idx])
            else:  # Return all three tensor components
                eps.extend(eps_[3 * i : 3 * i + 3])

        return eps, success

    def get_current_strains(self, gauge_ids):
        """
        This method computes the strain tensor at gauges from current state.

        Parameters
        ----------
        gauge_ids : list of int
            Array with gauge identification numbers

        Returns
        -------
        list of float
            Strain components
        bool
            Always True, unless the calculation failed
        """
        return self.compute_strains_from_displ(None, gauge_ids)

    def compute_rel_dist_from_displ(self, disp, ids):
        """
        This method computes the relative distance at sensors for the
        given displacement field.

        Parameters
        ----------
        disp : list of float
            Displacement field
        ids : list of int
            Array with function identification numbers

        Returns
        -------
        list of float
            Displacement components
        bool
            Always True, unless the calculation failed
        """

        n_dof_, disp_ = self._convert_c_double_array(disp)
        n_ids_, ids_ = self._convert_c_int_array(ids)

        n_ids = len(ids)
        rel_dis_ = (c_double * n_ids)()
        success = self._solver.getRelDisp(disp_, ids_, rel_dis_, n_dof_, n_ids_)

        rel_dis = [0.0] * n_ids
        rel_dis[:] = rel_dis_

        return rel_dis, success

    def compute_spring_var_from_displ(self, disp, ids):
        """
        This method computes one of the spring variables (length,
        deflection, force) at sensors for the given displacement field.

        Parameters
        ----------
        disp : list of float
            Displacement field
        ids : list of int
            Array with function identification numbers

        Returns
        -------
        list of float
            Displacement components
        bool
            Always True, unless the calculation failed
        """

        n_dof_, disp_ = self._convert_c_double_array(disp)
        n_ids_, ids_ = self._convert_c_int_array(ids)

        n_ids = len(ids)
        spr_var_ = (c_double * n_ids)()
        success = self._solver.getRespVars(disp_, ids_, spr_var_, n_dof_, n_ids_)

        spr_var = [0.0] * n_ids
        spr_var[:] = spr_var_

        return spr_var, success

    def compute_int_forces_from_displ(self, disp, ids):
        """
        This method computes beam section forces in triads for the
        given displacement field.

        There are 2 cases:

        | InternalForce component defined by an array
        |   dofi contains dof number in a certain direction i
        |   ids = [ beamID1, triadID1, dof1, beamID2, triadID2, dof2, ..... ]

        | InternalSectionForces defined by an array
        |   dofi contains always -1, 6 components are requested
        |   ids = [ beamID1, triadID1, -1, beamID2, triadID2, -1, ..... ]

        Parameters
        ----------
        disp : list of float
            Displacement field
        ids : list of int
            Array with beam and triad identification numbers

        Returns
        -------
        list of float
            Force components
        bool
            Always True, unless the calculation failed
        """

        n_dof_, disp_ = self._convert_c_double_array(disp)
        n_ids_, ids_ = self._convert_c_int_array(ids)

        n_items = n_ids_.value // 3
        n_items_ = c_int(n_items)
        forces_ = (c_double * (6 * n_items))()
        success = self._solver.getBeamForcesFromDisp(
            disp_, ids_, forces_, n_dof_, n_items_
        )

        iforce = 0
        forces = []
        for i in range(2, len(ids), 3):
            if ids[i] < 0:  # all 6 force components are requested
                forces.extend(forces_[iforce : iforce + 6])
            else:  # Force component for local dof ids[i] is requested
                forces.append(forces_[iforce + ids[i]])
            iforce += 6

        return forces, success

    def get_joint_spring_stiffness(self, bid):
        """
        Get joint spring stiffness coefficient(s).

        Parameters
        ----------
        bid : int
            Joint identification number (base ID)

        Returns
        -------
        list of float
            Spring stiffness coefficients for the joint DOFs
        bool
            Always True, unless the extraction failed
        """
        max_dof = 6  # max no. of joint dofs
        stiffc_ = (c_double * max_dof)()
        success = self._solver.getJointSprCoeff(stiffc_, self._convert_c_int(bid))

        stiffness = [0.0] * max_dof
        stiffness[:] = stiffc_

        return stiffness, success
