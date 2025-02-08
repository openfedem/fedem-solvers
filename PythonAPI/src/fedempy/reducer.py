# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Python wrapper for the native Fedem FE part reducer.
Used for convenience in order to hide native type convertions.
"""

from ctypes import c_bool, c_char_p, c_int, cdll


class FedemReducer:
    """
    This class mirrors the functionality of the Fedem FE part reducer library
    (libfedem_reducer_core.so on Linux, fedem_reducer_core.dll on Windows).

    Parameters
    ----------
    lib_path : str
        Absolute path to the reducer shared object library
    solver_options : list of str, default=None
        List of command-line arguments passed to the reducer

    Methods
    -------
    solver_init:
        Initializes the command-line handler of the reducer
    run:
        Invokes the reducer
    """

    def __init__(self, lib_path, solver_options=None):
        """
        Constructor.
        Optionally initialize the reducer itself if solver_options is given.
        """
        # load the reducer library
        self._solver = cdll.LoadLibrary(lib_path)

        # set up return type for functions in the reducer library
        self._solver.initSolverArgs.restype = c_int
        self._solver.reducePart.restype = c_int

        # initialize the fedem reducer
        self.__first_part = c_bool(True)
        self.solver_init(solver_options)

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
            argv[:] = [arg[i].encode("utf-8") for i in range(argc)]
            return c_int(argc), argv
        if (
            hasattr(arg, "_type_")
            and hasattr(arg, "_length_")
            and isinstance(arg, c_char_p)
        ):
            return c_int(len(arg)), arg

        raise TypeError(f"Expected {list} or {c_char_p}, got {type(arg)}.")

    def solver_init(self, options):
        """
        This function initializes the command-line handler of the reducer.

        See the Fedem R8.0 Users Guide, Appendix C.2 for a complete list of all
        command-line arguments that may be specified and their default values.

        Parameters
        ----------
        options : list of str
            List of command-line arguments passed to the reducer

        Returns
        -------
        int
            Zero on success, a negative value indicates some error condition
        """
        if options is None:
            return 0  # do nothing if no solver options specified

        # The first option is expected to contain the path of the executable,
        # so insert a dummy name here if the first option starts with '-'
        if isinstance(options, list) and options[0][0] == "-":
            options.insert(0, "fedem_reducer")  # arbitrary name (not used)
        argc_, argv_ = self._convert_c_char_array(options)

        return self._solver.initSolverArgs(argc_, argv_, self.__first_part)

    def run(self, options=None):
        """
        This function invokes the reducer.

        Parameters
        ----------
        options : list of str, default=None
            List of command-line arguments passed to the reducer

        Returns
        -------
        int
            Zero on success, a negative value indicates some error condition
        """
        self.solver_init(options)
        status = self._solver.reducePart(self.__first_part, c_bool(False))
        self.__first_part = c_bool(False)

        return status
