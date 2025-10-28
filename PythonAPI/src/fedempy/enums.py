# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
This module contains some enum definitions for the ``fedempy`` package.
They are python equivalents to corresponding enum definitions in the C++ code
of the `FedemDB` shared object library (see the `fedem_mdb` repository).
The numerical value of each enum value in this module should therefore
not be changed without a similar change in the C++ equivalent,
to preserve consistency.
"""

from enum import Enum


class FmType(Enum):
    """
    This enumerator identifies the various object types a Fedem mechanism model
    may consist of. They are mapped onto a corresponding type ID value of the
    C++ classes in the `FedemDB.C` source file (see the `fedem_mdb` repository).
    """

    ALL = 0
    TRIAD = 1
    BEAM = 2
    FEPART = 3
    BEAM_PROP = 4
    MAT_PROP = 5
    JOINT = 6
    RIGID_JOINT = 7
    REVOLUTE_JOINT = 8
    BALL_JOINT = 9
    FREE_JOINT = 10
    PRISMATIC_JOINT = 11
    CYLINDRIC_JOINT = 12
    CAM_JOINT = 13
    LOAD = 14
    FUNCTION = 15
    SENSOR = 16
    AXIAL_SPRING = 17
    AXIAL_DAMPER = 18
    STRAIN_ROSETTE = 19
    USERDEF_ELEMENT = 20
    REFERENCE_PLANE = 21


class FmLoadType(Enum):
    """
    This enumerator identifies the two external load types a Fedem model.
    The values corresponds to the enum values in the C++ class ``FmLoad``.
    (see the `fedem_mdb` repository).
    """

    FORCE = 0
    TORQUE = 1


class FmDofStat(Enum):
    """
    This enumerator identifies the available DOF constraint types.
    The values corresponds to the ``FmHasDOFsBase::DOFStatus`` enum values
    (see the `fedem_mdb` repository).
    """

    FREE = 0
    FIXED = 1
    PRESCRIBED = 2
    FREE_DYN = 3
    SPRING = 4
    SPRING_DYN = 5


class FmDof(Enum):
    """
    This enumerator identifies the local DOF components in an object.
    The values corresponds to the ``FmIsMeasuredBase::SensorDof`` enum values
    (see the `fedem_mdb` repository).
    """

    TX = 0
    TY = 1
    TZ = 2
    RX = 3
    RY = 4
    RZ = 5
    LENGTH = 6
    MAX_PR = 13
    MIN_PR = 14
    SA_MAX = 15
    VMISES = 16
    GAGE_1 = 17
    GAGE_2 = 18
    GAGE_3 = 19


class FmVar(Enum):
    """
    This enumerator identifies the result quantities that may be measured.
    The values corresponds to the ``FmIsMeasuredBase::SensorEntity`` enum values
    (see the `fedem_mdb` repository).
    """

    POS = 0
    LOCAL_VEL = 1
    GLOBAL_VEL = 2
    LOCAL_ACC = 3
    GLOBAL_ACC = 4
    DISTANCE = 5
    VEL = 6
    ACC = 7
    REL_POS = 8
    SPR_ANG = 9
    SPR_DEFL = 10
    SPR_FORCE = 11
    DAMP_ANG = 12
    DAMP_VEL = 13
    DAMP_FORCE = 14
    LENGTH = 15
    DEFLECTION = 16
    FORCE = 17
    LOCAL_FORCE = 18
    GLOBAL_FORCE = 19
    STRAIN = 24
    STRESS = 25
