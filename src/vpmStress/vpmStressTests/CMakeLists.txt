# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( LIB_ID vpmStressTests )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( FTENV_WARNINGS AND LINUX AND NOT USE_INTEL_FORTRAN )
  string ( APPEND CMAKE_Fortran_FLAGS " -Wno-argument-mismatch" )
endif ( FTENV_WARNINGS AND LINUX AND NOT USE_INTEL_FORTRAN )

include_directories ( ${CMAKE_CURRENT_BINARY_DIR}/.. )

add_library ( ${LIB_ID} ../elStressModule.f90
                        ../strainAndStressUtils.f90
                        ../../vpmCommon/samModule.f90
                        recKindModule.f90 ffl.f90 )
target_link_libraries ( ${LIB_ID} FemLib )

if ( pFUnit_FOUND )
  # Add some unit tests using pFUnit (executed via ctest)
  enable_fortran_tests ()
  add_fortran_test ( test_thickshell testThickShell.pf ${LIB_ID} )
endif ( pFUnit_FOUND )
