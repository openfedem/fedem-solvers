# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( LIB_ID vpmLinAlgTests )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

# Avoid that the generated module files are moved to the common folder
unset ( CMAKE_Fortran_MODULE_DIRECTORY )

# Need to use FFlLib classes in the FTK namespace
string ( APPEND CMAKE_CXX_FLAGS " -DFF_NAMESPACE=FTK" )

add_library ( ${LIB_ID} testModels.C )
add_library ( ${LIB_ID}_F90 generate_system.f90
                        ../../../vpmReducer/samReducerModule.f90
                        ../../../vpmReducer/rigidModule.f90
                        ../../../vpmReducer/wavgmModule.f90 )
target_link_libraries ( ${LIB_ID} FFlLib )
target_link_libraries ( ${LIB_ID}_F90 FFlLib_F ${LIB_ID} vpmLinAlg )

if ( pFUnit_FOUND )
  # Add some unit tests using pFUnit (executed via ctest)
  enable_fortran_tests ()
  add_fortran_test ( test_sysmat testSysMat.pf vpmLinAlg )
  add_fortran_test ( test_sprshift testSprShift.pf ${LIB_ID} ${LIB_ID}_F90 )
endif ( pFUnit_FOUND )
