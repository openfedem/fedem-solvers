# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Python wrapper for the native fedem dynamics solver
which operates directly on fedem mechanism model files (`*.fmm`).

This module relies on the following environment variables:

| *FEDEM_MDB* = Full path to the Fedem model database shared object library
| *FEDEM_SOLVER* = Full path to the Fedem dynamics solver shared object library
| *FEDEM_REDUCER* = Full path to the Fedem reducer shared object library
| *VIS_EXPORTER* = Full path to the Fedem VTFx exporter shared object library
The third variable needs to be set only if FE model reduction is to be performed.
The fourth variable needs to be set only if a VTFx file is to be exported.
The first two variables are mandatory.

This module can also be launched directly, to run a specified model,
using the syntax

| ``python -m fedempy.fmm_solver -f mymodel.fmm [--reduce-fem] [--save-model] [-v mymodel.vtfx]``

This will then invoke the method :meth:`fmm_solver.FmmSolver.solve_all`
on the specified model-file (`mymodel.fmm`). If you already have a directory
populated with Fedem solver input files, you can run the solver directly on
those files by launching this module without arguments, i.e.:

| ``python -m fedempy.fmm_solver``
"""

from argparse import ArgumentParser
from ctypes import c_double
from os import environ, getcwd, path

from fedempy.exporter import Exporter
from fedempy.fmm import FedemModel, FmType
from fedempy.inverse import InverseSolver
from fedempy.reducer import FedemReducer
from fedempy.solver import FedemException, FedemSolver


def _solver_options(solver, rdbdir="", time_start=None, base_id=None):
    """
    Sets up the standard solver command-line options,
    taking into account the *.fco, *.fop and *.fao files, if they exist.
    """
    opts = ["-cwd", rdbdir, "-terminal", "-1"] if rdbdir else []
    for ext in ("fco", "fop", "fao"):
        option_file = solver + "." + ext
        if path.isfile(rdbdir + option_file):
            opts += ["-" + ext, option_file]

    if time_start:  # Add dynamics solver start time
        opts.append("-timeStart=" + str(time_start))

    if base_id:  # Add flag for in-core von Mises stress array
        opts.append("-partVMStress=2")

    return opts


def _count_lines(file_name):
    """
    Counts the number of lines in a file.
    """
    count = 0
    if path.isfile(file_name):
        with open(file_name, "r") as file:
            for _ in file:
                count += 1

    return count


def _get_resfile(rdbdir, solver_name):
    """
    Extracts the actual res-file name from the fop-file for a solver process.
    """
    fop_file = rdbdir + solver_name + ".fop"
    if path.isfile(fop_file):
        with open(fop_file, "r") as fopf:
            for line in fopf:
                if line.find("-resfile") == 0:
                    return rdbdir + line[10:-2]  # Remove the enclosing ""

    return rdbdir + solver_name + ".res"


def _print_res(file_name):
    """
    Prints the contents of the specified res-file, the first 100 lines
    and the last part of it from the first Error message occuring, if any.
    """
    print("     Check", file_name, "for error messages.")
    print("     Here is the file content:")
    with open(file_name, "r") as resf:
        print_line = 1
        for count, line in enumerate(resf):
            if line.find("Error :") == 0:
                print_line = 2
            elif count > 100 and print_line == 1:
                print_line = 0
                print("     . . .")
            if print_line > 0:
                print("%8d %s" % (count + 1, line.rstrip()))
    print("#### End of file", file_name, flush=True)


class FmmSolver(FedemSolver):
    """
    This subclass of FedemSolver adds the possibility
    to start a simulation on a specified fedem model file.

    Parameters
    ----------
    model_file : str, default=None
        Absolute path to the fedem model file to start the simulation on
    use_internal_state : bool, default=False
        If True, internal state arrays are allocated (for microbatching)
    keep_old_res : bool, default=False
        Option to not overwrite any existing res-file in the RDB directory

    Methods
    -------
    start:
        Starts a simulation on the specified model file
    close_model:
        Closes the currently open model
    solve_all:
        Starts and runs through a simulation on the specified model file
    set_input:
        Assigns sensor value to a tagged external function
    """

    def __init__(self, model_file=None, use_internal_state=False, keep_old_res=False):
        """
        Constructor.
        Optionally starts the simulation, if a model file is specified.
        """
        if "FEDEM_SOLVER" not in environ:
            raise FedemException("Environment variable FEDEM_SOLVER not defined")
        if "FEDEM_MDB" not in environ:
            raise FedemException("Environment variable FEDEM_MDB not defined")

        super().__init__(environ["FEDEM_SOLVER"], None, use_internal_state)
        self._model = FedemModel(environ["FEDEM_MDB"])
        self._reducer = None
        self._vtfx = None
        self._func_map = {}
        self._c_transf = None
        self._c_deform = {}
        self._c_stress = {}
        self._prev_time = 0
        self._time_incr = 0

        # Start the Fedem dynamics solver
        status = self.start(model_file, keep_old_res)
        if status < 0:
            raise FedemException(
                "Failed to start Fedem on " + model_file + f" ({status})"
            )

    def _init_func_map(self):
        """
        Initializes the tag -> channel_index map for external functions.
        """
        chn = 1
        tag = self._model.fm_get_func_tag(chn)
        while tag is not None:
            if tag:
                self._func_map[tag] = chn
            chn += 1
            tag = self._model.fm_get_func_tag(chn)

    def _reduce_fe_model(self, do_reduce):
        """
        Reduces all FE parts in the currently open model,
        which are not already redused.

        Parameters
        ----------
        do_reduce : bool
            Skip the FE part reduction if False

        Returns
        -------
        int
            Number of FE parts that was reduced, negative on error
        """

        num_reduced = 0
        if not do_reduce:
            return num_reduced

        if self._reducer is None:
            if "FEDEM_REDUCER" not in environ:
                return num_reduced
            self._reducer = FedemReducer(environ["FEDEM_REDUCER"])

        # Get list (of base Id) of the FE parts in the model
        fe_parts = self._model.fm_get_objects(FmType.FEPART)
        for base_id in fe_parts:
            # Check if this part is reduced and create reducer input files if not
            rdbdir = self._model.fm_write_reducer(base_id)
            if len(rdbdir) > 0:
                # Run the FE part reducer
                num_reduced += 1
                if self._reducer.run(_solver_options("fedem_reducer", rdbdir)) != 0:
                    print(" *** Reduction failure for FE part", base_id)
                    _print_res(_get_resfile(rdbdir, "fedem_reducer"))
                    return -num_reduced

                print("   * FE Part", base_id, "successfully reduced", flush=True)
                self._model.fm_sync_part(base_id)

        if num_reduced > 0:
            print("#### FE model reduction done.", flush=True)
        return num_reduced

    def _open_vtfx_exporter(self, vtfx_file, model_file):
        """
        Opens the VTFx file exporter and write the FE models to it.

        Parameters
        ----------
        vtfx_file : str
            Absolute path of the VTFx output file
        model_file : str
            Fedem model file name to derive the case name from

        Returns
        -------
        list of int
            Base ID of FE parts for which stress recovery will be performed
        """

        if vtfx_file is None:
            return None

        if "VIS_EXPORTER" not in environ:
            print("  ** Environment variable VIS_EXPORTER not defined, no VTFx export.")
            return None

        # Get list (of base Id) of the FE parts in the model
        # and the associated FE data files
        base_ids = []
        fem_parts = []
        vis_parts = []
        all_parts = self._model.fm_get_objects(FmType.FEPART)
        for base_id in all_parts:
            fem_file, recover = self._model.fm_get_femodel(base_id)
            if recover == -1:  # Generic part with visualization file only
                vis_parts.append(
                    {
                        "path": fem_file,
                        "name": path.splitext(path.basename(fem_file))[0],
                        "base_id": base_id,
                    }
                )
            elif recover >= 0:  # This is an FE part
                fem_parts.append(
                    {
                        "path": fem_file,
                        "name": path.splitext(path.basename(fem_file))[0],
                        "base_id": base_id,
                        "recovery": recover > 0,
                    }
                )
            if recover > 0:  # Stress recovery will be performed for this FE part
                base_ids.append(base_id)

        # Create the VTFx exporter object
        self._vtfx = Exporter(
            fem_parts,
            vis_parts,
            environ["VIS_EXPORTER"],
            vtfx_file,
            path.splitext(path.basename(model_file))[0],
        )

        return base_ids

    def _init_result_buffers(self, base_ids):
        """
        Allocates state arrays to pass results to the VTFx exporter module.

        Parameters
        ----------
        base_ids : list of int
            Base ID of FE parts for which stress recovery will be performed

        """
        self._c_transf = (c_double * self.get_transformation_state_size())()
        self._c_deform = {}
        self._c_stress = {}
        for bid in base_ids:
            def_size = self.get_part_deformation_state_size(bid)
            if def_size > 0:
                str_size = self.get_part_stress_state_size(bid)
                self._c_deform[bid] = (c_double * def_size)()
                self._c_stress[bid] = (c_double * (str_size if str_size > 0 else 0))()

    def start(
        self,
        model_file,
        keep_old_res=False,
        close_model=False,
        reduce_fem=False,
        state_data=None,
        gauge_data=None,
        extf_input=None,
        time_start=None,
        vtfx_file=None,
    ):
        """
        Starts a simulation on the specified model file.

        Parameters
        ----------
        model_file : str
            Absolute path of the Fedem model file to run the solver on
        keep_old_res : bool, default=False
            Option to not overwrite any existing res-file in the RDB directory
        close_model : bool, default=False
            If True, release the model from memory before solver start
        reduce_fem : bool, default=False
            If True, and the model contains FE parts, they will be reduced
            before the solver is started, unless the reduced matrix files
            already exist
        state_data : list of float, default=None
            Complete state vector to restart simulation from
        gauge_data : list of float, default=None
            Initial strain gauge values for restart
        extf_input : list of float, default=None
            Initial external function values, for initial equilibrium iterations
        time_start : float, default=None
            Optional start time of simulation, override setting in model file
        vtfx_file : str, default=None
            Absolute path of the VTFx output file for visualization in GLview

        Returns
        -------
        int
            Zero on success, otherwise negative
        """
        if model_file is None:
            return 0  # do nothing if no model file specified

        # Count the number of lines in the existing log-file
        log_file = model_file.replace(".fmm", ".log")
        this_log = _count_lines(log_file)

        # Local function for printing the model log-file when shit happens.
        def error_exit(status):
            self._model.fm_close(True)
            print("     Here is the model log:")
            with open(log_file, "r") as logf:
                for count, line in enumerate(logf):
                    if count >= this_log:
                        print(line.rstrip())
            print("#### End of log for", model_file, flush=True)
            return status

        # Read the fedem model file into the internal datastructure
        if not self._model.fm_open(model_file):
            print(" *** Failed to open model file", model_file)
            print("     in working directory", getcwd())
            return error_exit(-97)

        print("   * Model file", model_file, "successfully opened")
        print("     Number of Triads:", self._model.fm_count(FmType.TRIAD))
        print("     Number of Joints:", self._model.fm_count(FmType.JOINT))
        print("     Number of Beams:", self._model.fm_count(FmType.BEAM))
        print("     Number of Parts:", self._model.fm_count(FmType.FEPART))
        print("     Total number of mechanism objects:", self._model.fm_count())

        # Check for FE model reduction
        num_red = self._reduce_fe_model(reduce_fem)
        if num_red < 0:
            return error_exit(-98)

        # Initialize the VTFx file exporter
        base_ids = self._open_vtfx_exporter(vtfx_file, model_file)

        # Create the results database file structure
        # populated with the solver input files
        rdbdir = self._model.fm_write_solver(keep_old_res)
        if len(rdbdir) < 1:
            print(" *** Failed to write solver input.")
            return error_exit(-99)

        # Initialize the tag -> channel_index map for external functions
        self._init_func_map()

        # Optionally close the model before starting the solver.
        # The model file can then not be updated with the new results.
        if close_model:
            self.close_model(num_red > 0)

        # Start the dynamics solver
        options = _solver_options("fedem_solver", rdbdir, time_start, base_ids)
        status = self.solver_init(options, None, state_data, gauge_data, extf_input)
        if status < 0:
            if not close_model:
                self._model.fm_close(True)
            print(" *** Failed to start the dynamics solver.")
            _print_res(_get_resfile(rdbdir, "fedem_solver"))
        elif base_ids is not None:
            self._init_result_buffers(base_ids)

        return status

    def _export_step(self):
        """
        Exports results for current step to the VTFx file, if any.
        """
        if not self._vtfx or self.have_results() == 0:
            return

        self.save_transformation_state(self._c_transf)
        for bid, deformation in self._c_deform.items():
            self.save_part_state(bid, deformation, self._c_stress[bid])

        # Find the (largest) time increment, for setting the animation frame rate
        ctime = self.get_current_time()
        delta = ctime - self._prev_time
        self._prev_time = ctime
        if delta > self._time_incr:
            self._time_incr = delta

        self._vtfx.do_step(ctime, self._c_transf, self._c_deform, self._c_stress)

    def close_model(self, save=False, remove_singletons=False, fringe_max=-1):
        """
        Closes the currently open model.

        Parameters
        ----------
        save : bool, default=False
            If True, the model file is updated based on current result files
        remove_singletons : bool, default=False
            If True, heap-allocated singelton objects are also released
        fringe_max : float, default=-1
            Max fringe range value for VTFx output.
            If less than zero, the range will be set automatically.

        Returns
        -------
        bool
            True on success, otherwise False
        """
        if save:
            print("   * Saving updated model")
            status = self._model.fm_save()
        else:
            status = True

        self._model.fm_close(not status or remove_singletons)

        if self._vtfx:
            self._vtfx.clean(self._time_incr, fmax=fringe_max)

        return status

    def solve_all(
        self,
        model_file,
        save_model=True,
        reduce_fem=False,
        vtfx_file=None,
        fringe_max=-1,
    ):
        """
        Starts and runs through a simulation on the specified model file.

        Parameters
        ----------
        model_file : str
            Absolute path of the Fedem model file to run the solver on
        save_model : bool, default=True
            If True, save the model file with new results when solver finished
        reduce_fem : bool, default=False
            If True, and the model contains FE parts, they will be reduced
            before the solver is started, unless the reduced matrix files
            already exist
        vtfx_file : str, default=None
            Absolute path of the VTFx output file for visualization in GLview
        fringe_max : float, default=-1
            Max fringe range value for VTFx output.
            If less than zero, the range will be set automatically.

        Returns
        -------
        int
            Zero on success, otherwise negative
        """

        if model_file is None:
            # Run the solver directly in current working directory,
            # if it contains some solver input files (fco, fop, etc).
            # Otherwise, do nothing.
            options = _solver_options("fedem_solver")
            return self.run_all(options) if options else 0

        print("\n#### Running dynamics solver on", model_file)
        ierr = self.start(model_file, True, False, reduce_fem, vtfx_file=vtfx_file)
        if ierr < 0:
            print(f" *** Solver failed to start ({ierr}).")
            return ierr

        # Run through the entire time series
        while self.solve_next():
            self._export_step()

        if self.solver_done() == 0 and self.ierr.value == 0:
            print("     Time step loop OK, solver closed")
            self.close_model(save_model, False, fringe_max)
        else:
            self.close_model(False, True, fringe_max)
            print(" *** Dynamics solver failed", self.ierr.value)

        return self.ierr.value

    def set_input(self, func_tag, value=None):
        """
        This method may be used prior to the solve_next call, to assign a
        sensor value from a physical twin to the specified actuator or load in
        the model, identified by the argument func_tag.

        Parameters
        ----------
        func_tag : str
            Tag of the external function to be assigned new value
        value : float, default=None
            The value to be assigned

        Returns
        -------
        bool
            Always True, unless the specified function is not found
        """
        if func_tag in self._func_map:
            return self.set_ext_func(self._func_map[func_tag], value)
        return False


class FmmInverse(FmmSolver, InverseSolver):
    """
    This class augments FmmSolver with inverse solution capabilities.

    Parameters
    ----------
    model_file : str
        Absolute path to the fedem model file to start the simulation on
    config : dictionary
        Inverse solver configuration
    use_internal_state : bool, default=False
        If True, internal state arrays are allocated (for microbatching)
    keep_res : bool, default=False
        Option to not overwrite any existing res-file in the RDB directory
    """

    def __init__(self, model_file, config, use_internal_state=False, keep_res=False):
        """
        Constructor.
        Optionally starts the simulation, if a model file is specified.
        """
        # Initialize the dynamics solver
        FmmSolver.__init__(self, model_file, use_internal_state, keep_res)

        # Initialize the inverse solver
        InverseSolver.__init__(self, self, config)


if __name__ == "__main__":
    parser = ArgumentParser(description="Python Fedem Solver")
    parser.add_argument(
        "-f", "--model-file", required=True, help="Fedem model file to open"
    )
    parser.add_argument(
        "-s", "--save-model", action="store_true", help="Save updated model file"
    )
    parser.add_argument(
        "-r", "--reduce-fem", action="store_true", help="Reduce model before solve"
    )
    parser.add_argument(
        "-v", "--vtfx-file", help="Ceetron VTFx output file for visualization"
    )
    parser.add_argument(
        "-m", "--fringe-max", type=float, default=-1, help="Max von Mises fringe value"
    )
    model = FmmSolver()
    model.solve_all(**vars(parser.parse_args()))
