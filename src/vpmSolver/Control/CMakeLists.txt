# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID Control )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )
if ( BUILD_CONTROL_AS_DLL )
  set ( LIB_ID fedem_solver_ctrl )
endif ( BUILD_CONTROL_AS_DLL )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( USE_INTEL_FORTRAN )
  set ( BLA_VENDOR Intel10_64lp )
endif ( USE_INTEL_FORTRAN )
find_package ( LAPACK REQUIRED )

file ( GLOB F77_SOURCE_FILES *.f )
file ( GLOB F90_SOURCE_FILES *.f90 )


# Build and install

if ( BUILD_CONTROL_AS_DLL )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_SHARED" )
  add_library ( ${LIB_ID} SHARED ${F77_SOURCE_FILES} ${F90_SOURCE_FILES} )
  target_link_libraries ( ${LIB_ID} LINK_PRIVATE ${LAPACK_LIBRARIES} )
  install ( TARGETS ${LIB_ID}
            RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
            LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" )
else ( BUILD_CONTROL_AS_DLL )
  add_library ( ${LIB_ID} STATIC ${F77_SOURCE_FILES} ${F90_SOURCE_FILES} )
  target_link_libraries ( ${LIB_ID} vpmUtilities_F90 ${LAPACK_LIBRARIES} )
endif ( BUILD_CONTROL_AS_DLL )
