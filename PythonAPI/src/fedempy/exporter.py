# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Python wrapper for the native VTFX exporter.
"""

from ctypes import byref, c_bool, c_char_p, c_double, c_float, c_int, cdll
from os import path, remove
from subprocess import run


class ExporterException(Exception):
    """
    General exception-type for exporter exceptions.

    Parameters
    ----------
    rc : int, default=0
        Return code
    """

    def __init__(self, rc=0):
        """
        Constructor.
        """
        super().__init__({"Failure": f"Return code {rc}"})


class Exporter:
    """
    This class provides functionality for exporting fedem animations
    to Ceetron CUG format.

    Parameters
    ----------
    fe_parts : dict
        Dictionary of finite element parts to include in visualization
    vis_parts : dict
        Dictionary of visualization/vrml parts to include in visualization
    lib_path : str
        Absolute path to the shared object library vtfxExporter.so
    vtfx_path : str, default=None
        Absolute path to vtfx-file, use a temporary file if None
    case_name : str, default="Case"
        Case identifier

    Methods
    -------
    do_step:
        Executes a step of visualization export with the provided input data
    clean:
        Converts temporary generated vtfx-file to CUG database and cleans up
    """

    def __init__(self, fe_parts, vis_parts, lib_path, vtfx_path=None, case_name="Case"):
        """
        Constructor.
        Initializes the internal datastructure of the shared object library,
        and loads the finite element and visualization parts into memory.
        """
        self._initialize(lib_path, vtfx_path, None, case_name)
        self._fem_parts = []
        self._vis_parts = []

        # Add FE parts
        for fe_part in fe_parts:
            if path.isfile(fe_part["path"]):
                self._add_fe_part(
                    str.encode(fe_part["path"]),
                    str.encode(fe_part["name"]),
                    fe_part["base_id"],
                    fe_part.get("surface_only", True),
                    fe_part.get("include_spiders", False),
                )
                if fe_part["recovery"]:
                    self._fem_parts.append(fe_part["base_id"])
                else:
                    self._vis_parts.append(fe_part["base_id"])

        # Add visualization parts
        for vis_part in vis_parts:
            if path.isfile(vis_part["path"]):
                self._add_visualization_part(
                    str.encode(vis_part["path"]),
                    str.encode(vis_part["name"]),
                    vis_part["base_id"],
                )
                self._vis_parts.append(vis_part["base_id"])

    def _get_part_object(self, part_id):
        elements = self._get_elements(part_id)
        types = self._get_element_types(part_id)
        nodes = self._get_nodes(part_id)
        t_matrix = self._get_transformation_matrix(part_id)
        return {
            "nodes": nodes,
            "edges": elements,
            "elementTypes": types,
            "transformation": t_matrix,
        }

    def _export_geometry(self):
        parts = []
        for part_id in self._fem_parts:
            parts.append(self._get_part_object(part_id))
        for part_id in self._vis_parts:
            parts.append(self._get_part_object(part_id))

        return {"parts": parts}

    def do_step(self, time, transformation_in, deformation_in, stress_in):
        """
        Execute a step of visualization export with the provided input data.

        Parameters
        ----------
        time : float
            Simulation time
        transformation_in : list of c_double
            Array of transformation data from the fedem solver
        deformation_in : dict
            Arrays of part deformation data from the fedem solver
        stress_in : dict
            Arrays of part stress data from the fedem solver
        """
        self._set_transformation_data(transformation_in)
        for base_id in self._fem_parts:
            self._set_part_deformation_data(base_id, deformation_in[base_id])
            self._set_part_stress_data(base_id, stress_in[base_id])

        rc = self.lib_exporter.doStep(c_double(time))

        if rc != 0:
            raise ExporterException(rc)

    def clean(self, time_step, lib_dir=None, out_dir=None, fmax=-1.0):
        """
        Clears all data about included FE-parts and frs-files,
        and removes all transformation, stress and deformation results.
        Closes the vtfx-file and converts it to a CUG database.

        Parameters
        ----------
        time_step : float
            Time step to use for animation speed
        lib_dir : str, default None
            Absolute path for folder where to write CUG log-file
            If None, current working directory will be used
        out_dir : str, default=None
            Absolute path for output folder of the CUG database
            If None, CUG data is not generated and the vtfx-file is retained
        fmax : float, default=-1
            Set fringe range to [0,fmax] if fmax is larger than 0.
            If less than zero, the range will be set automatically.
        """
        print("   * Writing VTFX-file", self.vtfx_path, flush=True)
        rc = self.lib_exporter.finish(c_double(time_step), c_double(0), c_double(fmax))
        if rc == 0 and out_dir and path.isfile(self.cug_path):
            logfile = lib_dir + "/Cug.log" if path.isdir(lib_dir) else "./Cug.log"
            with open(logfile, "w") as outfile:
                print("   * Creating CUG database in", out_dir, flush=True)
                rc = run(
                    [self.cug_path, self.vtfx_path, out_dir], check=True, stdout=outfile
                ).returncode

        if rc != 0:
            print(f" *** Failed ({rc})")
        elif out_dir:
            remove(self.vtfx_path)

    def _initialize(self, lib_path, vtfx_path, cug_path, case_name):
        """
        Initialization
        """
        self.lib_exporter = cdll.LoadLibrary(lib_path)
        self.lib_exporter.initialize(c_bool(False))
        self.vtfx_path = "/tmp/temp.vtfx" if vtfx_path is None else vtfx_path
        self.cug_path = "/usr/bin/CugComposer" if cug_path is None else cug_path
        rc = self.lib_exporter.setVtfxPath(
            c_char_p(self.vtfx_path.encode("utf-8")),
            c_char_p(case_name.encode("utf-8")),
        )
        if rc != 0:
            raise ExporterException(rc)

    def _add_fe_part(
        self, ftl_path, part_name, base_id, surface_only=True, include_spiders=False
    ):
        """
        Adds the FE-Part defined in the FEDEM ftl-file at ftl_path,
        with the specified part_name, and set up remapping of elements
        to the formats required by the Ceetron vtfx format.

        Parameters
        ----------
        ftl_path : str
            Absolute path to the ftl-file containing the node
            and element definitions for the FE-part.
        part_name : str
            Name given to the FE-part in the exported vtfx-file.
        base_id : int
            The base_id of the part in the FEDEM model.
        surface_only : bool, default=True
            Set to True if all geometry should be converted to
            surface quad and triangle elements.
        include_spiders : bool, default=False
            Set to True to also include spider elements and beams
        """
        if not path.isfile(ftl_path):
            raise ExporterException(1)

        rc = self.lib_exporter.addFEPart(
            c_char_p(ftl_path),
            c_char_p(part_name),
            c_int(base_id),
            c_bool(surface_only),
            c_bool(include_spiders),
            c_float(1.0),
            c_float(1.0),
            None,
        )
        if rc != 0:
            raise ExporterException(rc)

    def _add_visualization_part(self, file_path, part_name, base_id):
        """
        Adds the visualization part defined in the specified vrml-file,
        with the specified part_name, and set up remapping of elements
        to the formats required by the Ceetron vtfx format.

        Parameters
        ----------
        file_path : str
            Absolute path to the visualization file containing the node
            and element definitions.
        part_name : str
            Name given to the FE-part in the exported vtfx-file.
        base_id : int
            The base_id of the part in the FEDEM model.
        """
        if not path.isfile(file_path):
            raise ExporterException(1)

        rc = self.lib_exporter.addVisualizationPart(
            c_char_p(file_path),
            c_char_p(part_name),
            c_int(base_id),
            c_float(1.0),
            c_float(1.0),
            None,
        )
        if rc != 0:
            raise ExporterException(rc)

    def _get_number_of_elements(self, base_id):
        """
        Returns the total number of finite elements for the FE-part in the ftl-file.

        Parameters
        ----------
        baseId : int
            The FEDEM BaseID of the FE-part.

        Returns
        -------
        int
            The number of elements.
        """
        return self.lib_exporter.getNumberOfElements(c_int(base_id))

    def _get_number_of_element_nodes(self, base_id):
        """
        Returns the total number of element-nodes for the FE-part
        with the specified baseId.

        Parameters
        ----------
        baseId : int
            The FEDEM BaseID of the FE-part.

        Returns
        -------
        int
            The number of element nodes.
        """
        return self.lib_exporter.getNumberOfElementNodes(c_int(base_id))

    def _get_number_of_nodes(self, base_id):
        """
        Returns the total number of unique nodes for the FE-part
         with the specified baseId.

        Parameters
        ----------
        baseId : int
            The FEDEM BaseID of the FE-part.

        Returns
        -------
        int
            The number of nodes.
        """
        return self.lib_exporter.getNumberOfNodes(c_int(base_id))

    def _get_element_types(self, base_id, element_types=None):
        """
        Fills the array given as input with the element types for all elements
        of the FE-part with the specified baseId.
        Each type is an integer-value corresponding to the number of nodes
        in the element(3 for triangles, 4 for quads)

        Parameters
        ----------
        baseId : int
            The FEDEM BaseId of the FE-Part.
        elementTypes : list of c_int, default=None
            The array of integers that will be filled with the element-types.
            If none, a list will be allocated inside the function.
            Size of list must be equal to value from _get_number_of_elements().

        Returns
        -------
        list of ints
            List of element types
        """
        if not element_types:
            element_types = (c_int * self._get_number_of_elements(base_id))()
        rc = self.lib_exporter.getElementTypes(c_int(base_id), byref(element_types))

        if rc != 0:
            raise ExporterException(rc)
        return element_types[0:]

    def _get_elements(self, base_id, elements=None):
        """
        Fills the array given as input with the node-indices of all
        elements of the FE-part with the specified baseId.

        Parameters
        ----------
        baseId : int
            The FEDEM BaseID of the FE-part.
        elements : list of c_int, default=None
            List of integers that will be filled with the element-node indices.
            If None, a list will be allocated inside the function.
            The size must be equal to the sum of all nodes for the different elements.

        Returns
        -------
        List of ints
            List of node indices for the elements
        """
        if not elements:
            elements = (c_int * self._get_number_of_element_nodes(base_id))()
        rc = self.lib_exporter.getElements(
            c_int(base_id), byref(elements), c_bool(False)
        )

        if rc != 0:
            raise ExporterException(rc)
        return elements[0:]

    def _get_nodes(self, base_id, nodes=None):
        """
        Fills the array given as input with the x,y and z values
        for all the nodes of the FE-part with the specified baseId.

        Parameters
        ----------
        baseId : int
            The FEDEM BaseID of the FE-part.
        nodes : list of c_double, default=None
            List of doubles that will be filled with the node coordinates.
            If None, a list will be allocated inside the function.
            The size must be equal to three times the value from _get_number_of_nodes().

        Returns
        -------
        List of doubles
            List of node coordinates
        """
        if not nodes:
            nodes = (c_double * (self._get_number_of_nodes(base_id) * 3))()
        rc = self.lib_exporter.getNodes(c_int(base_id), byref(nodes))

        if rc != 0:
            raise ExporterException(rc)
        return nodes[0:]

    def _get_transformation_matrix(self, base_id, transformation=None):
        """
        Fills the array given as input with the 16 components of the
        transformation matrix of the FE-part with the
        specified baseId, for the current time-step.
        Row dominant. Final row is translation

        Parameters
        ----------
        baseId : int
            The FEDEM BaseID of the part to get results from.
        transformation : list of c_double, default=None
            List of doubles that will be filled with the transformation matrix.
            If none, a list will be allocated inside the function.
            Size of array must be equal to 16.

        Returns
        -------
        List of doubles
            List containing the transformation matrix
        """
        if not transformation:
            transformation = (c_double * 16)()
        rc = self.lib_exporter.getTransformationMatrix(
            c_int(base_id), byref(transformation)
        )

        if rc != 0:
            raise ExporterException(rc)
        return transformation[0:]

    def _get_deformation_vector(self, base_id, deformation=None):
        """
        Fills the array given as input with the updated x,y and z values in the
        current time-step, for all the nodes of
        the FE-part with the specified baseId.

        Parameters
        ----------
        baseId : int
            The FEDEM BaseID of the part to get results from.
        deformation : List of c_double, default=None
            List of doubles that will be filled with
            x,y,z coordinates of all nodes of the deformed part.
            If None, a list will be allocated inside the function.
            The size must be equal to three times the value from _get_number_of_nodes().

        Returns
        -------
        List of doubles
            List containing the nodal displacements
        """
        if not deformation:
            deformation = (c_double * (self._get_number_of_nodes(base_id) * 3))()
        rc = self.lib_exporter.getDeformationVector(c_int(base_id), byref(deformation))

        if rc != 0:
            raise ExporterException(rc)
        return deformation[0:]

    def _get_stress_vector(self, base_id, stress=None):
        """
        Fills the array given as input with the von Mises stress values in the
        current time-step, for all the nodes of the FE-part in the ftl-file.

        Parameters
        ----------
        baseId : int
            The FEDEM BaseID of the part to get results from.
        stress : List of c_double, default=None
            List of doubles that will be filled with
            von Mises stress value of all nodes of the deformed part.
            If none, a list will be allocated inside the function.
            The size must be equal to the value from _get_number_of_element_nodes().

        Returns
        -------
        List of doubles
            List containing the stress values for each node
        """
        if not stress:
            stress = (c_double * self._get_number_of_nodes(base_id))()
        rc = self.lib_exporter.getStressVector(c_int(base_id), byref(stress))

        if rc != 0:
            raise ExporterException(rc)
        return stress[0:]

    def _set_transformation_data(self, data):
        sz = len(data)
        rc = self.lib_exporter.setTransformationData(c_int(sz), data)

        if rc != 0:
            raise ExporterException(rc)

    def _set_part_deformation_data(self, base_id, data):
        sz = len(data)
        rc = self.lib_exporter.setPartDeformationData(c_int(base_id), c_int(sz), data)

        if rc != 0:
            raise ExporterException(rc)

    def _set_part_stress_data(self, base_id, data):
        sz = len(data)
        rc = self.lib_exporter.setPartStressData(c_int(base_id), c_int(sz), data)

        if rc != 0:
            raise ExporterException(rc)
