# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Python implementation of inverse solution methods with Fedem.
"""

from copy import deepcopy
from os import environ, path

from numpy import array, c_, delete, dot, sqrt, transpose, vstack, zeros
from numpy.linalg import inv, lstsq, solve

from fedempy.enums import FmType
from fedempy.log_conf import get_logger
from fedempy.solver import FedemException, FedemSolver

try:
    from scipy import linalg

    have_sci_py = True
except ImportError:
    have_sci_py = False

# log file at the execution place
logger = get_logger("fedemRun.log")


class InverseException(FedemException):
    """
    General exception type for inverse solver exceptions.
    Used to generalize error messages from the FedemSolver methods.

    Parameters
    ----------
    method_name : str
        Name of the method that detected an error
    ierr : int, default=None
        Error flag value that is embedded into the error message
    """

    def __init__(self, method_name, ierr=None):
        """
        Constructor. Forwards to the parent class constructor.
        """
        error_message = "FedemSolver." + method_name + "() failure"
        if ierr is None or ierr == 0:
            super().__init__(error_message)
        else:
            super().__init__(
                error_message + f" ({ierr}). Check fedem_solver.res for details."
            )


class InverseSolver:
    """
    This class handles the inverse solution through proper methods.
    It accesses the Fedem model through the provided FedemSolver instance.

    Parameters
    ----------
    solver : FedemSolver
        The Fedem dynamics solver instance
    config : dictionary
        Inverse solver configuration

    Methods
    -------
    run_inverse_dyn:
        Performs the inverse solution (dynamic case)
    run_inverse_fedem:
        Performs the inverse static solution using the internal inverse solver
    run_inverse:
        Performs the inverse solution (static case)
    """

    def __init__(self, solver, config):
        """
        Constructor. Initializes the object.
        """
        # create dicts
        self.internal_equations = config.get("internal_equations", {})
        if config is not None and self.internal_equations:
            logger.info("initialisation list started")
            self._init_equations(solver)
            logger.info("initialisation list finished\n")

        self.solver = solver

        # Memory for dynamic inverse solution
        self.u_vec = None  # displacement vector
        self.ud_vec = None  # velocity vector
        self.udd_vec = None  # acceleration vector
        self.fa_vec = None  # force vector for the HHT implementation

        self.internal_force_mat = None  # initial force matrix
        self.loop_nr = 0  # loop number over simulation

    def _init_equations(self, solver):  # NOSONAR
        """
        Find internal equation number related to triad_id and dof
        """

        def __get_equations(solver, obj_id):
            """
            Return the equation numbers associated with given object(s).
            """
            if not isinstance(obj_id, str):
                return solver.get_equations(obj_id)
            meqn = []
            base_ids = solver._model.fm_get_objects(tag=obj_id)
            for bid in base_ids:
                meqn.extend(solver.get_equations(bid))
            return meqn

        dof_set = {"tx": 0, "ty": 1, "tz": 2, "rx": 3, "ry": 4, "rz": 5}
        dof_set_rev = {"rz": 0, "tz": 1}

        # allowed set of equations
        eq_list = [
            "unknown_fm",
            "unknown_f",
            "known_x",
            "known_intF",
            "known_secF",
            "known_relD",
            "known_eps",
            "known_sprD",
            "known_sprF",
            "known_Fx",
        ]

        # if 'baseID' in dict, change it to 'triadID'
        # ignore key unknown_fm
        ignore_key = [eq_list[0]]
        for k, eq_items in self.internal_equations.items():
            if k not in ignore_key:
                for item in eq_items:
                    if "baseID" in item.keys():
                        item["triadID"] = item.pop("baseID")

        # defined set of equations (collects the keys in the list)
        # and use the same order  as mentioned in the yaml file
        self.eq_list_def = [*self.internal_equations]

        # declare empty lists
        self.modes = []
        self.int_force_list = []
        self.sec_force_list = []
        self.rel_dist_list = []
        self.strain_tensor_list = []
        self.defl_spring_list = []
        self.frc_spring_list = []

        # Manage eigenvalue calculation, only once or every time step
        self.use_initial_eigen_vec = False

        # Select which eigensolver to use.
        # The default (0) is to use Fedem's internal Lanczos solver.
        # Notice that the other options involve expanding the system matrices
        # into full dense matrices, which will require higher memory usage and
        # longer computation time. Therefore, use for smaller systems only.
        self.modes_solver = 0  # 1: DSYGVX, 2: DGGEVX, 3: scipy.linalg.eigh

        for item in eq_list:
            if item in self.internal_equations:
                if item == eq_list[0]:
                    if "initial" in self.internal_equations[item]:
                        self.use_initial_eigen_vec = True
                        self.modes = self.internal_equations[item]["initial"]
                        logger.info("Eigenmodes/eigenvectors are calculated only once")
                    elif "update" in self.internal_equations[item]:
                        self.modes = self.internal_equations[item]["update"]
                    else:
                        self.modes = self.internal_equations[item]
                    logger.info("Using modes: %s" % self.modes)
                    print("Mode numbers: ", self.modes)
                    if "solver" in self.internal_equations[item]:
                        self.modes_solver = [
                            "LANCZOS",
                            "DSYGVX",
                            "DGGEVX",
                            "SCIPY_EIGH",
                        ].index(self.internal_equations[item]["solver"])
                        logger.info("Using eigensolver: %s" % self.modes_solver)
                elif item in (eq_list[1], eq_list[2], eq_list[9]):
                    # branch for "unknown_f", "known_x", "known_Fx"
                    item_len = len(self.internal_equations[item])
                    for i in range(item_len):
                        if "triadID" in self.internal_equations[item][i]:
                            con_id = self.internal_equations[item][i]["triadID"]
                            con_dof_set = dof_set
                        else:
                            con_id = self.internal_equations[item][i]["revJID"]
                            con_dof_set = dof_set_rev

                        try:
                            val = int(self.internal_equations[item][i]["dof"])
                            v_list = []
                            while val:
                                digit = val % 10
                                v_list.append(digit - 1)
                                # remove last digit from number (as integer)
                                val //= 10
                            eq_num = __get_equations(solver, con_id)
                            for idx, v_item in enumerate(v_list):
                                v_list[idx] = eq_num[v_item]

                            # append dict in reverse order
                            self.internal_equations[item][i]["eqNum"] = v_list[::-1]
                        except ValueError as value_error:
                            dof = con_dof_set[self.internal_equations[item][i]["dof"]]
                            eq_num = __get_equations(solver, con_id)
                            if len(eq_num) == 0:
                                raise ValueError(
                                    "Check if triadID is correct, "
                                    + "dof is not found or may be fixed"
                                ) from value_error
                            if eq_num[dof] < 0:
                                raise ValueError(
                                    "Check if triadID is correct, "
                                    + "dof is likely dependent"
                                ) from value_error

                            eq_num = eq_num[dof]

                            self.internal_equations[item][i]["eqNum"] = eq_num
                elif item == eq_list[3]:
                    self.int_force_list = self.internal_equations[item]
                elif item == eq_list[4]:
                    self.sec_force_list = self.internal_equations[item]
                elif item == eq_list[5]:  # known_relD - relative distance
                    self.rel_dist_list = self.internal_equations[item]
                elif item == eq_list[6]:  # known_eps - plain strain tensor
                    self.strain_tensor_list = self.internal_equations[item]
                elif item == eq_list[7]:  # known_sprD - spring deflection
                    self.defl_spring_list = self.internal_equations[item]
                elif item == eq_list[8]:  # known_sprF - spring force
                    self.frc_spring_list = self.internal_equations[item]

    def _inverse_get_boundary_conditions(self, one_based=False):
        """
        Extracts boundary conditions related to measurements and forces.
        """

        def __get_equation_numbers(equation_set):
            """
            Extracts equation numbers from given dictionary.
            """
            eq_list = []
            for ieq in equation_set:
                eq_num = ieq["eqNum"]
                if isinstance(eq_num, list):
                    eq_list.extend(eq_num)
                else:
                    eq_list.append(eq_num)
            if not one_based:  # fortran to python convention
                eq_list = [x - 1 for x in eq_list]
            return eq_list

        x_def = None
        g_def = None

        item = self.internal_equations.get("known_x", None)
        item = self.internal_equations.get("known_Fx", item)
        if item:
            x_def = __get_equation_numbers(item)
            logger.info("--> content of x_def: %s" % x_def)

        item = self.internal_equations.get("unknown_f", None)
        if item:
            g_def = __get_equation_numbers(item)
            logger.info("--> content of g_def: %s" % g_def)

        return x_def, g_def

    def convert_rev_joint_force(self, data):
        """
        Modify the data set for spring force input for revolute joints
        with defined spring forces
        Conversation from force to lenght x=F/k (const. stiffness assumption)

        Parameters
        ----------
        data : list of float
            Input function/data values

        Returns
        -------
        list of int
            revolute joint ID's and their modified input data
        """
        dof_set_rev = {"rz": 0, "tz": 1}

        rev_joints = []

        item = self.internal_equations.get("known_Fx")
        if item:
            for idx, val in enumerate(item):
                jid = val.get("revJID", -1)
                dof = val.get("dof", "rz")
                if jid > 0 and dof in ("rz", "tz"):
                    dof = dof_set_rev[dof]
                    rev_joints.append(jid)
                    stiff = self.solver.get_joint_spring_stiffness(jid)[0][dof]
                    print("spring stiffness: ", stiff)
                    data[idx] /= stiff

        return rev_joints

    @staticmethod
    def _inverse_build_mat_from_ndef(k_mat, mat_def):
        """
        Building position matrices for active dofs (sensor/force)
        Indicating related dof by 1
        """
        mat = zeros((len(k_mat), len(mat_def)))
        for idx, val in enumerate(mat_def):
            mat[val, idx] = 1

        return mat

    def _init_mem(self, ndof, incs=2):
        """
        Memory initialisation for dynamic inverse solution
        """
        self.u_vec = zeros((ndof, incs), float)  # displacement vector
        self.ud_vec = zeros((ndof, incs), float)  # velocity vector
        self.udd_vec = zeros((ndof, incs), float)  # acceleration vector
        self.fa_vec = zeros((ndof, incs), float)  # RHS force vector

    @staticmethod
    def _newmark_coefficients(h, dt):
        """
        Parameter calculation for Newmark and HHT algorithm
        """
        alpha = 0.25 * (1.0 - h) * (1.0 - h)
        beta = 0.5 * (1.0 - (2.0 * h))

        a0 = 1.0 / (alpha * dt * dt)
        a1 = beta / (alpha * dt)
        a2 = 1.0 / (alpha * dt)
        a3 = (1.0 / (2.0 * alpha)) - 1.0
        a4 = (beta / alpha) - 1.0
        a5 = (dt / 2.0) * ((beta / alpha) - 2.0)

        return a0, a1, a2, a3, a4, a5

    def run_inverse_dyn(self, inp_data, out_def):
        """
        Inverse solution driver (dynamic case).

        Parameters
        ----------
        inp_data : list of float
            Input function values
        out_def : list of int
            User Ids of the functions to evaluate the response for

        Returns
        -------
        list of float
            Evaluated response variables
        """
        if self.loop_nr == 0:
            self._init_mem(self.solver.get_system_size())

        t_0 = self.solver.get_current_time()

        # run start step
        do_continue = self.solver.start_step()
        if self.solver.ierr.value < 0:  # Simulation failure
            raise InverseException("start_step", self.solver.ierr.value)
        if not do_continue:  # Reached the end of simulation
            return None

        t1 = self.solver.get_current_time()

        q_vec, ok = self.solver.get_external_force_vector()
        if not ok:
            raise InverseException("get_external_force_vector")

        x_def, g_def = self._inverse_get_boundary_conditions()
        print("x_def and g_def: ", x_def, "  ", g_def)

        # simulation parameters
        h = -0.1  # default Newmark alpha value
        dt = t1 - t_0
        a0, a1, a2, a3, a4, a5 = self._newmark_coefficients(h, dt)

        k_mat = self.solver.get_stiffness_matrix()[0]
        m_mat = self.solver.get_mass_matrix()[0]
        c_mat = self.solver.get_damping_matrix()[0]
        n_mat = self.solver.get_newton_matrix()[0]
        z_mat = self._inverse_build_mat_from_ndef(n_mat, x_def)

        # assign measurement data/sensor data
        measurements = inp_data

        if self.loop_nr > 0:
            # define storage indices
            i0 = (self.loop_nr - 1) % 2
            i1 = (self.loop_nr) % 2

            # displacements v0, velocities v1 and accelerations v2
            # are currently calulated 'online'
            v1 = (
                a1 * self.u_vec[:, i0]
                + a4 * self.ud_vec[:, i0]
                + a5 * self.udd_vec[:, i0]
            )
            v2 = (
                a0 * self.u_vec[:, i0]
                + a2 * self.ud_vec[:, i0]
                + a3 * self.udd_vec[:, i0]
            )

            cv = dot(c_mat, v1)
            ma = dot(m_mat, v2)

            fac = self.fa_vec[:, i0]

            # equation: kh*u = F + cv + ma (where kh = n_mat)
            # 1) calculate the dynamic part from the former solution step
            # 2) measurements (sensor values) are split into an external force part
            #    (generalized forces) and dynamic part (cv and ma are calculated
            #    from the step before)
            # 3) reduce the measurements by the dynamic part
            # 4) calculate the external force vector (generalized forces) by
            #    calling _inverse_core()
            # 5) generalized forces are scaled by parameter h (fedem alignment)
            csc = (1.0 + h) * (cv) + ma - h * fac
            uc = solve(n_mat, csc)
            uc0 = dot(transpose(z_mat), uc)

            # call inverse_core method,
            # dynamic part is subtracted from the measurements
            fsc = self._inverse_core(n_mat, x_def, g_def, measurements - uc0, q_vec)[1]

            # force (for fedem input)
            f_pos = (1.0 / (1.0 + h)) * fsc

            # force for displacement, velocity and acceleration calculation
            f_dyn = fsc + csc

            # new displacements
            u_n = solve(n_mat, f_dyn)

            # update accelerations (u_ddn) and velocities (u_dn)
            u_ddn = a0 * u_n - v2
            u_dn = a1 * u_n - v1

            fan = f_pos - dot(c_mat, u_dn) - dot(k_mat, u_n)

            # store displacements, velocities and accelerations
            self.u_vec[:, i1] = u_n
            self.ud_vec[:, i1] = u_dn
            self.udd_vec[:, i1] = u_ddn
            self.fa_vec[:, i1] = fan

        else:
            # first increment uses static solution
            f_pos = self._inverse_core(k_mat, x_def, g_def, measurements, q_vec)[1]

        # subtract constant force vector
        f_pos -= q_vec

        # update right hand side vector
        if not self.solver.add_rhs_vector(f_pos):
            raise InverseException("add_rhs_vector")

        # equilibrium iterations (fedem)
        self.solver.finish_step()
        if self.solver.ierr.value != 0:
            raise InverseException("finish_step", self.solver.ierr.value)

        # increase time step loop number
        self.loop_nr += 1

        # return output function values
        return self.solver.get_functions(out_def)

    def _inverse_core(self, k_mat, x_def, g_def, x_vec, f0=None):
        """
        Central inverse algorithm
        """
        b_mat = self._inverse_build_mat_from_ndef(k_mat, x_def)
        f_mat = self._inverse_build_mat_from_ndef(k_mat, g_def)

        if f0 is None:
            f0 = zeros(len(k_mat))  # constant forces

        inv_k = inv(k_mat)
        bt_inv_k = dot(b_mat.transpose(), inv_k)
        c_val = dot(bt_inv_k, f_mat)
        b0 = dot(bt_inv_k, f0)

        opt_alpha = lstsq(c_val, (x_vec - b0), rcond=-1)[0]
        print("scaling: ", opt_alpha)
        force_eval = dot(f_mat, opt_alpha) + f0
        pos_eval = dot(inv_k, force_eval)  # displacement vector

        return pos_eval, force_eval

    def run_inverse_fedem(self, inp_data, out_def):
        """
        This method uses fedem's inverse solution in fortran

        Parameters
        ----------
        inp_data : list of float
            Input function values
        out_def : list of int
            User Ids of the functions to evaluate the response for

        Returns
        -------
        list of float
            Evaluated response variables
        """
        x_def, g_def = self._inverse_get_boundary_conditions(True)
        out, _ = self.solver.solve_inverse(inp_data, x_def, g_def, out_def)
        if self.solver.ierr.value < 0:
            raise InverseException("solve_inverse", self.solver.ierr.value)

        return out

    def _build_sensor_ids(self, sensor_type):  # NOSONAR
        """
        Store sensor IDs
        """

        def __get_base_id(obj_id, obj_type):
            """
            Converts an object tag into corresponsing base Id.
            """
            if isinstance(obj_id, str):
                base_ids = self.solver._model.fm_get_objects(obj_type, obj_id)
                if len(base_ids) == 1:
                    obj_id = base_ids[0]
                elif len(base_ids) > 1:
                    print(f"Multiple objects of type {obj_type} with tag {obj_id}")
                    obj_id = 0
                else:
                    print(f"No objects of type {obj_type} with tag {obj_id}")
                    obj_id = 0

            if obj_id > 0:
                return obj_id

            raise ValueError(
                "Check if sensor_id is correct, "
                + "sensor_id is not found or may be not defined"
            )

        sensor_id = []
        nr = 0

        if (sensor_type == "relative") and (self.rel_dist_list):
            item_type = FmType.TRIAD
            sensor_list = self.rel_dist_list
            identifier = "relID"
        elif (sensor_type == "strain") and (self.strain_tensor_list):
            item_type = FmType.STRAIN_ROSETTE
            sensor_list = self.strain_tensor_list
            identifier = "epsID"
        elif (sensor_type == "force") and (self.int_force_list):
            item_type = FmType.BEAM
            sensor_list = self.int_force_list
            identifier = "beamID"
        elif (sensor_type == "section") and (self.sec_force_list):
            item_type = FmType.BEAM
            sensor_list = self.sec_force_list
            identifier = "beamID"
        elif (sensor_type == "deflection") and (self.defl_spring_list):
            item_type = FmType.JOINT
            sensor_list = self.defl_spring_list
            identifier = "deflID"
        elif (sensor_type == "springForce") and (self.frc_spring_list):
            item_type = FmType.JOINT
            sensor_list = self.frc_spring_list
            identifier = "frcID"
        else:
            print("no sensor type specified")
            return sensor_id, nr

        for item in sensor_list:
            item_len = len(item)
            sensor_id.append(__get_base_id(item[identifier], item_type))

            if identifier == "beamID":
                sensor_id.append(__get_base_id(item["triadID"], FmType.TRIAD))
                if item_len >= 3:  # case internal force
                    try:  # check if d_set an int or a string
                        d_set = int(item["dof"])
                        v_list = []
                        while d_set:
                            digit = d_set % 10
                            # switch from fortran to python notation
                            v_list.append(digit - 1)
                            d_set //= 10
                        # insert last element
                        sensor_id.append(v_list[-1])
                        nr += 1
                        # loop backwards, v_list contains digits in reverse order
                        for k in range(len(v_list) - 2, -1, -1):
                            sensor_id.append(sensor_id[-3])  # append beamID again
                            sensor_id.append(sensor_id[-3])  # append triadID again
                            sensor_id.append(v_list[k])  # append dof
                            nr += 1
                    except ValueError:
                        dof_set = {"tx": 0, "ty": 1, "tz": 2, "rx": 3, "ry": 4, "rz": 5}
                        dof = dof_set[item["dof"]]
                        sensor_id.append(dof)
                        nr += 1
                else:  # case section forces
                    sensor_id.append(-1)
                    nr += 1
            elif identifier == "epsID":
                if item_len >= 2:  # strain tensor component is defined
                    try:  # check if e_set an int or a string
                        e_set = int(item["strain"])
                        v_list = []
                        while e_set:
                            digit = e_set % 10
                            # switch from fortran to python notation
                            v_list.append(digit - 1)
                            e_set //= 10
                        # insert last element
                        sensor_id.append(v_list[-1])
                        nr += 1
                        # loop backwards, v_list contains digits in reverse order
                        for k in range(len(v_list) - 2, -1, -1):
                            sensor_id.append(sensor_id[-2])  # append epsID again
                            sensor_id.append(v_list[k])  # append direction
                            nr += 1
                    except ValueError:
                        eps_set = {"ex": 0, "ey": 1, "exy": 2}
                        eps = eps_set[item["strain"]]
                        sensor_id.append(eps)
                        nr += 1
                else:
                    sensor_id.append(-1)
                    nr += 1
            elif identifier in ("relID", "deflID", "frcID"):
                nr += 1

        print(f"No. of sensors for {sensor_type}: {nr}, sensor_id={sensor_id}")

        return sensor_id, nr

    @staticmethod
    def _inverse_disp_sensor(dim, x_def, u_vec, lhs, rhs, pos):
        """
        Standard inverse solution for displacement sensor input
        """
        b_mat = zeros((dim, len(x_def)))

        for idx, val in enumerate(x_def):
            b_mat[val, idx] = 1

        c_val = dot(b_mat.transpose(), u_vec)

        # remove non-scalable part from the measurements
        for k in range(len(x_def)):
            rhs[k] -= c_val[k][-1]

        lhs = vstack([lhs, c_val[:, :-1]])
        pos += len(x_def)

        return lhs, rhs, pos

    def _inverse_strain_sensor(self, u_vec, lhs, rhs, pos):
        """
        Inverse solution for strain sensor input (e.g. strain gage measurements)
        """
        gage_id, nr_gauges = self._build_sensor_ids("strain")

        n_3c = gage_id.count(-1)  # number of 3 component tensors
        n_1c = nr_gauges - n_3c  # number of 1 component tensors
        nr_c = 3 * n_3c + n_1c

        # no. of generalized load cases (without gravity load displacement vector)
        ng = len(u_vec[0]) - 1

        # gathering unit load strains into a matrix
        lh = zeros((nr_c, ng))

        # loop over generalized forces
        for j in range(ng):
            lh[:, j], ok = self.solver.compute_strains_from_displ(u_vec[:, j], gage_id)
            if not ok:
                raise InverseException("compute_strains_from_displ")

        # const force part (e.g. gravity)
        eps, ok = self.solver.compute_strains_from_displ(u_vec[:, -1], gage_id)
        if not ok:
            raise InverseException("compute_strains_from_displ")

        # update measurements based on strains, subtract constant part
        for k in range(nr_c):
            rhs[pos + k] -= eps[k]

        lhs = vstack([lhs, lh])
        pos += nr_c

        return lhs, rhs, pos

    def _inverse_section_forces_sensor(self, u_vec, lhs, rhs, pos):
        """
        Inverse solution for section forces (N,Qy,Qz,Mx,My,Mz)
        6 components solution
        """
        beam_ids, nr_sec = self._build_sensor_ids("section")

        # no. of generalized load cases (without gravity load displacement vector)
        ng = len(u_vec[0]) - 1

        # gathering sectional forces (from unit loads) into a matrix
        lh = zeros((nr_sec * 6, ng))

        # loop over generalized forces
        for j in range(ng):
            lh[:, j], ok = self.solver.compute_int_forces_from_displ(
                u_vec[:, j], beam_ids
            )
            if not ok:
                raise InverseException("compute_int_forces_from_displ")

        c_fg, ok = self.solver.compute_int_forces_from_displ(u_vec[:, -1], beam_ids)
        if not ok:
            raise InverseException("compute_int_forces_from_displ")

        # remove constant load from unit force matrix
        # update measurements/calculations based on internal forces
        for j in range(nr_sec * 6):
            rhs[pos + j] -= c_fg[j]

        lhs = vstack([lhs, lh])
        pos += nr_sec * 6

        return lhs, rhs, pos

    def _inverse_rel_dist_sensor(self, u_vec, lhs, rhs, pos):
        """
        Inverse solution for relative distance change
        """
        eng_ids, _ = self._build_sensor_ids("relative")

        # no. of generalized load cases (without gravity load displacement vector)
        ng = len(u_vec[0]) - 1

        # define size of the lhs matrix
        lh = zeros((len(eng_ids), ng))

        # loop over generalized forces (no. of columns)
        for j in range(ng):
            lh[:, j], ok = self.solver.compute_rel_dist_from_displ(u_vec[:, j], eng_ids)
            if not ok:
                raise InverseException("compute_rel_dist_from_displ")

        # relative displacements from constant loads (e.g. gravity)
        r_dg, ok = self.solver.compute_rel_dist_from_displ(u_vec[:, -1], eng_ids)
        if not ok:
            raise InverseException("compute_rel_dist_from_displ")

        for k in range(len(eng_ids)):
            rhs[pos + k] -= r_dg[k]

        lhs = vstack([lhs, lh])
        pos += len(eng_ids)

        return lhs, rhs, pos

    def _inverse_spring_var(self, u_vec, lhs, rhs, pos, spr_var):
        """
        Inverse solution for spring variables
        """
        ids, _ = self._build_sensor_ids(spr_var)

        # no. of generalized load cases (without gravity load displacement vector)
        ng = len(u_vec[0]) - 1

        # define size of the lhs matrix
        lh = zeros((len(ids), ng))

        # loop over generalized forces (no. of columns)
        for j in range(ng):
            lh[:, j], ok = self.solver.compute_spring_var_from_displ(u_vec[:, j], ids)
            if not ok:
                raise InverseException("compute_spring_var_from_displ")

        # spring forces from constant loads (e.g. gravity)
        r_dg, ok = self.solver.compute_spring_var_from_displ(u_vec[:, -1], ids)
        if not ok:
            raise InverseException("compute_spring_var_from_displ")

        for k in range(len(ids)):
            rhs[pos + k] -= r_dg[k]

        lhs = vstack([lhs, lh])
        pos += len(ids)

        return lhs, rhs, pos

    def _inverse_int_force_sensor(self, u_vec, lhs, rhs, pos):
        """
        Inverse solution for known internal force
        """
        ids, nr_b = self._build_sensor_ids("force")

        # no. of generalized load cases (without gravity load displacement vector)
        ng = len(u_vec[0]) - 1

        # gathering forces (from unit loads) into a matrix
        lh = zeros((nr_b, ng))

        # loop over generalized forces
        for j in range(ng):
            lh[:, j], ok = self.solver.compute_int_forces_from_displ(u_vec[:, j], ids)
            if not ok:
                raise InverseException("compute_int_forces_from_displ")

        # beam forces from const forces (e.g. gravity)
        c_fg, ok = self.solver.compute_int_forces_from_displ(u_vec[:, -1], ids)
        if not ok:
            raise InverseException("compute_int_forces_from_displ")

        # remove constant load from unit force matrix
        # update measurements/calculations based on internal forces
        for j in range(nr_b):
            rhs[pos + j] -= c_fg[j]

        lhs = vstack([lhs, lh])
        pos += nr_b

        return lhs, rhs, pos

    @staticmethod
    def _unit_load(dim, g_def):
        """
        Force vector based on unit loads (generalized forces)
        """
        # Prepare matrix for generalized forces
        f_mat = zeros((dim, len(g_def)))

        for idx, val in enumerate(g_def):
            f_mat[val, idx] = 1.0

        return f_mat

    @staticmethod
    def _mode_load_scipy(solver, modes):
        """
        Force vector based on natural frequency shapes
        """
        k_mat, ok = solver.get_stiffness_matrix()
        if not ok:
            raise InverseException("get_stiffness_matrix")

        m_mat, ok = solver.get_mass_matrix()
        if not ok:
            raise InverseException("get_mass_matrix")

        # check is matrix m_mat positive definite (via Cholesky factorization)
        print("Check Mass Matrix for positive definiteness:")
        try:
            linalg.cholesky(m_mat)
        except linalg.LinAlgError:
            print("Mass matrix is not positive definite - check input")
            logger.info(
                "Mass matrix will be modified, because matrix is not positive definite"
            )
            for i in range(m_mat.shape[0]):
                if m_mat[i, i] < 1.0e-15:
                    m_mat[i, i] = 1.0e-15

        # take into account modes up to max mode number in the array modes
        up_to_mode = max(modes) - 1
        print("calculate modes up to mode No.: ", up_to_mode)
        logger.info("Calculate modes up to mode No.: %s" % up_to_mode)

        # calculate natural frequencies (e) and natural frequency modes (u)
        (e_val, e_vec) = linalg.eigh(k_mat, m_mat, eigvals=(0, up_to_mode))
        print("Eigenvalues (low, high):", e_val[0], e_val[-1])
        logger.info("Eigenvalues: %s" % e_val)

        dim = solver.get_system_size()
        n_m = len(modes)
        f_vec = zeros((dim, n_m))
        fm = zeros((up_to_mode + 1, n_m))
        for idx, imode in enumerate(modes):
            fm[imode - 1, idx] = 1.0

        # Transform force vector from modal space (fm) to nodal space (F)
        # via F = (E^-1)^T*fm
        # The eigenvector matrix is an orthogonal matrix,
        # where transpose and inverse delivers the same result.
        # The above expression can therefore be simplified to F = E*fm

        f_vec = dot(e_vec, fm)
        logger.info("Modal force vector calculated")

        return f_vec

    @staticmethod
    def _mode_load(solver, modes, use_lapack):
        """
        Force vector based on natural frequency shapes.
        New and faster implementation, using Fedem's internal eigenvalue solver.
        """
        # take into account modes up to max mode number in the array modes
        n_modes = max(modes)
        print("Calculating the modes: ", modes)
        logger.info("Calculate the first %s eigenmodes." % n_modes)

        # calculate natural frequencies (e_val)
        # and the associated mode shapes (e_vec)
        e_val, e_vec, ok = solver.solve_modes(n_modes, False, use_lapack)
        if solver.ierr.value < 0 or not ok:
            raise InverseException("solve_modes", solver.ierr.value)
        if e_val is None:
            print(f"Unable for calculate eigenvalues at t={solver.get_current_time()}")
            print(f"Terminating dynamics solver ({solver.solver_done(print_res=True)})")
            print("Please check the fedem_solver.res file content above.")
            raise InverseException("solve_modes", solver.ierr.value)

        print("Eigenvalues (low, high):", e_val[0], e_val[-1])
        logger.info("Eigenvalues: %s" % e_val)

        # Pick eigenvectors as given by the modes indices
        n_dim = solver.get_system_size()
        n_mod = len(modes)
        f_vec = zeros((n_dim, n_mod))
        for idx, imode in enumerate(modes):
            f_vec[:, idx] = e_vec[imode - 1]
        logger.info("Modal force vector calculated")

        return f_vec

    def run_inverse(self, inp_data, out_def):  # NOSONAR
        """
        Collector for different inverse methods.

        Parameters
        ----------
        inp_data : list of float
            Input function values
        out_def : list of int
            User Ids of the functions to evaluate the response for

        Returns
        -------
        list of float
            Evaluated response variables
        """
        # run start step
        logger.info("================================")
        logger.info("Running step start %s" % self.loop_nr)
        do_continue = self.solver.start_step()
        if self.solver.ierr.value < 0:  # Simulation failure
            raise InverseException("start_step", self.solver.ierr.value)
        if not do_continue:  # Reached the end of simulation
            return None

        logger.info("Getting updated stiffness matrix")
        k_mat, ok = self.solver.get_stiffness_matrix()
        if not ok:
            raise InverseException("get_stiffness_matrix")

        # external force vector
        logger.info("Getting external force vector")
        q_vec, ok = self.solver.get_external_force_vector()
        if not ok:
            raise InverseException("get_external_force_vector")

        logger.info("Setting boundary conditions")
        x_def, g_def = self._inverse_get_boundary_conditions()
        print("x_def and g_def: ", x_def, "  ", g_def)

        # displacement field based on unit loads/modal forces
        # F matrix (generalized fore part) will not change during the simulation,
        # the modal force is calculated every step (default)
        # if use_initial_eigen_vec = False
        # otherwise only once at the beginning due performance reasons.
        f_mat = None
        n_dim = self.solver.get_system_size()

        if self.internal_force_mat is None or not self.use_initial_eigen_vec:
            if g_def is not None:
                f_mat = self._unit_load(n_dim, g_def)
                logger.info("Generalized force built")

            if self.modes:
                if self.modes_solver == 3 and have_sci_py:
                    f_vec = self._mode_load_scipy(self.solver, self.modes)
                else:
                    f_vec = self._mode_load(self.solver, self.modes, self.modes_solver)
                if f_mat is None:
                    f_mat = f_vec
                else:
                    f_mat = c_[f_mat, f_vec]
                logger.info("Modal force vector built")

            if self.internal_force_mat is None:
                self.internal_force_mat = deepcopy(f_mat)

        if self.use_initial_eigen_vec:
            logger.info("Store force vector (generalized and/or modal) at initial step")
            f_mat = deepcopy(self.internal_force_mat)

        # add external force vector to F (last column)
        f_mat = c_[f_mat, q_vec]

        # calculate displacement vector u by solving eq. k_mat*u=F
        logger.info("Calculating generalized displacements")
        u_vec = solve(k_mat, f_mat)

        # build rhs vector from measurements (copy)
        logger.info("Building RHS vector on sensor data")
        rhs = [0.0] * len(inp_data)
        rhs[:] = inp_data

        # position in rhs vector
        pos = 0

        # create a 1xg_def matrix, initialized with 0
        lhs = [[0] * (len(f_mat[0]) - 1) for _ in range(1)]

        # building equation system
        for eq_def in self.eq_list_def:
            if eq_def in ("known_x", "known_Fx"):
                lhs, rhs, pos = self._inverse_disp_sensor(
                    n_dim, x_def, u_vec, lhs, rhs, pos
                )
                logger.info(
                    "Displacement sensor, locations/pos [%s/%s]" % (len(lhs) - 1, pos)
                )
            elif eq_def == "known_eps":
                lhs, rhs, pos = self._inverse_strain_sensor(u_vec, lhs, rhs, pos)
                logger.info(
                    "Strain sensor, locations/pos [%s/%s]" % (len(lhs) - 1, pos)
                )
            elif eq_def == "known_secF":
                lhs, rhs, pos = self._inverse_section_forces_sensor(
                    u_vec, lhs, rhs, pos
                )
                logger.info(
                    "Section force sensor, locations/pos [%s/%s]" % (len(lhs) - 1, pos)
                )
            elif eq_def == "known_relD":
                lhs, rhs, pos = self._inverse_rel_dist_sensor(u_vec, lhs, rhs, pos)
                logger.info(
                    "Rel. dist sensor, locations/pos [%s/%s]" % (len(lhs) - 1, pos)
                )
            elif eq_def == "known_intF":
                lhs, rhs, pos = self._inverse_int_force_sensor(u_vec, lhs, rhs, pos)
                logger.info(
                    "Internal force sensor, locations/pos [%s/%s]" % (len(lhs) - 1, pos)
                )
            elif eq_def == "known_sprD":
                lhs, rhs, pos = self._inverse_spring_var(
                    u_vec, lhs, rhs, pos, "deflection"
                )
                logger.info(
                    "Spring deflection sensor, locations/pos [%s/%s]"
                    % (len(lhs) - 1, pos)
                )
            elif eq_def == "known_sprF":
                lhs, rhs, pos = self._inverse_spring_var(
                    u_vec, lhs, rhs, pos, "springForce"
                )
                logger.info(
                    "Spring force sensor, locations/pos [%s/%s]" % (len(lhs) - 1, pos)
                )

        # remove first row from matrix
        lhs = delete(lhs, (0), axis=0)

        # check number of equations and pointer position (pos)
        if (len(lhs) != len(rhs)) or (len(rhs) != pos):
            print("InCorrect dimension of equation system lhs: ", len(lhs))
            print("dimension of rhs: ", len(rhs))
            print("pos of position pointer: ", pos)
            logger.info(
                "Equation system has incorrect dimensions: [%s/%s]" % (len(lhs), pos)
            )

        print("Size of the equation system: ", len(lhs), "x", len(lhs[0]))

        # compute scaling factor and force vector
        alpha = lstsq(lhs, rhs, rcond=-1)[0]
        force = dot(f_mat[:, :-1], alpha)

        print("scaling factor: ", alpha)
        logger.info("Scaling factors: %s" % alpha)

        # update right hand side vector
        logger.info("Setting RHS for Fedem solver")
        if not self.solver.add_rhs_vector(force):
            raise InverseException("add_rhs_vector")

        print("force vector added")

        # equilibrium iterations (fedem)
        self.solver.finish_step()
        if self.solver.ierr.value != 0:
            raise InverseException("finish_step", self.solver.ierr.value)

        print("equlibrium iterations done")
        logger.info("Equlibrium iterations finished")

        # increment counter
        self.loop_nr += 1

        logger.info("-- Step/Cycle finished --\n")

        # return output function values
        return self.solver.get_functions(out_def)

    @staticmethod
    def _calc_eigen_values(solver):
        """
        Calulates eigenvalues and eigenvectors (test routine, only for testing)
        """
        if not have_sci_py:
            raise FedemException("FedemRun._calc_eigen_values() requires scipy")

        k_mat, ok = solver.get_stiffness_matrix()
        if not ok:
            raise InverseException("get_stiffness_matrix")

        m_mat, ok = solver.get_mass_matrix()
        if not ok:
            raise InverseException("get_mass_matrix")

        print("All eigenvalues/eigenvectors calculated")
        (eig_vals, eig_vecs) = linalg.eig(k_mat, m_mat)
        if any(eig_vals) < -1.0e-5:
            print("Warning: check k_mat, m_mat - negative eigenvalues")

        omega = array(sqrt(abs(eig_vals)) / 2.0 / 3.14159)
        order = omega.ravel().argsort()

        ndof = len(eig_vals)

        # put eigenvalues into correct order into fn
        fn = zeros(ndof)
        for i in range(0, ndof):
            fn[i] = omega[order[i]]

        # mass normalisation V_T*M*V
        dd = dot(eig_vecs.T, dot(m_mat, eig_vecs))
        for i in range(0, ndof):
            nf = sqrt(dd[i, i])
            for j in range(0, ndof):
                eig_vecs[j, i] /= nf

        # sort eigenvectors into MS
        ms = zeros((ndof, ndof))
        for i in range(0, ndof):
            ms[0:ndof, i] = eig_vecs[0:ndof, order[i]]

        print("Eigenvalues: ", fn)
        print("Eigenvectors: ", ms)


class FedemRun(FedemSolver, InverseSolver):
    """
    This class augments FedemSolver with inverse solution capabilities.

    Parameters
    ----------
    wrkdir : str
        Current working directory for the fedem dynamics solver
    config : dictionary
        Content of yaml input file
    """

    def __init__(self, wrkdir, config):
        """
        Constructor. Initializes the object.
        """
        if "FEDEM_SOLVER" not in environ:
            raise FedemException("Environment variable FEDEM_SOLVER not defined")

        # Set up the standard solver command-line options,
        # taking into account the *.fco, *.fop and *.fao files, if they exist.
        args = ["-cwd", wrkdir, "-terminal", "-1"]
        for ext in ("fco", "fop", "fao"):
            option_file = "fedem_solver." + ext
            if path.isfile(wrkdir + "/" + option_file):
                args += ["-" + ext, option_file]

        # Initialize the dynamics solver
        FedemSolver.__init__(self, environ["FEDEM_SOLVER"], args)

        # Initialize the inverse solver
        InverseSolver.__init__(self, self, config)
