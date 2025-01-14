# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( LIB_ID FemLibTest )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( pFUnit_FOUND )
# Add some unit tests using pFUnit (executed via ctest)
  enable_fortran_tests ()
  add_fortran_test ( test_ANDES testANDES.pf FemLib )
  add_fortran_test ( test_FqtSh testFqtS.pf  FemLib )
  add_fortran_test ( test_SCqtS testSCqtS.pf FemLib )
  add_fortran_test ( test_CongT testCongruenceTrans56.pf FemLib )
  add_fortran_test ( test_MPRO3 testMPRO30.pf FemLib )
endif ( pFUnit_FOUND )

add_executable ( ${LIB_ID} linkTest.f90 )
target_link_libraries ( ${LIB_ID} FemLib )
