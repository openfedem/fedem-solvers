# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Note: The tests here are put in a separate sub-folder,
# since they depend on a dummy environmenttypemodule
# (in file dummyModules.f90) which is different than the
# same module used by the other tests in the parent folder.
# This way a separate mod-file is created for this dummy module
# which does not disturb the other tests if run separately.

set ( LIB_ID moreTests )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( pFUnit_FOUND )
  enable_fortran_tests ()
  add_fortran_test ( test_ramp testSmoothRamp.pf vpmCommon_F90 FFaFunctionLib_F90 )
  target_sources ( test_ramp PUBLIC
                   ../../engineRoutinesModule.f90
                   ../../functionTypeModule.f90
                   ../../sensorTypeModule.f90
                   ../../springTypeModule.f90
                   ../../damperTypeModule.f90
                   dummyModules.f90 )
endif ( pFUnit_FOUND )
