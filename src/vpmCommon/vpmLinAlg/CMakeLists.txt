# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID vpmLinAlg )
set ( LIB_ID_LIST )
if ( BUILD_TESTS )
  set ( LIB_ID_LIST vpmLinAlgTests )
endif ( BUILD_TESTS )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

# Find external libraries for the 3rd-party linear equation solver

option ( USE_LINALG_GSF "Build with GSF linear equation solver" OFF )
mark_as_advanced ( USE_LINALG_GSF )
if ( USE_LINALG_GSF )
  find_package ( GSF )
  if ( GSF_LIBRARY )
    message ( STATUS "INFORMATION : Building ${LIB_ID} with GSF support" )
  endif ( GSF_LIBRARY )
endif ( USE_LINALG_GSF )

if ( USE_INTEL_FORTRAN )
  set ( BLA_VENDOR Intel10_64lp )
endif ( USE_INTEL_FORTRAN )
find_package ( LAPACK REQUIRED )
if ( LAPACK_LIBRARIES MATCHES mkl_intel )
  message ( STATUS "Using Intel Math Kernal Library" )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_MKL" )
endif ( LAPACK_LIBRARIES MATCHES mkl_intel )

if ( TARGET SPR_I8 )
  message ( STATUS "Using 64-bit integer version of the SPR library" )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_SPR_INT8" )
endif ( TARGET SPR_I8 )

# Subfolder handling

set ( SUBFOLDER_LIST ${LIB_ID_LIST} )
if ( NOT GSF_LIBRARY )
# When the 3rd-party equation solver is not available,
# add a dummy library containing the expected symbols.
# This library is compiled in a sub-folder, to avoid that
# its generated module files are mixed with the others.
  list ( INSERT SUBFOLDER_LIST 0 vpmLinAlgDummy )
endif ( NOT GSF_LIBRARY )
foreach ( FOLDER ${SUBFOLDER_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${SUBFOLDER_LIST} )


## Fortran90 source files
set ( F90_FILE_LIST asmExtensionModule
                    feLanczosModule
                    matExtensionModule
                    solExtensionModule
                    sprExtensionModule
                    sysMatrixTypeModule
    )
## Fortran90 callbacks from the linear equation solver
set ( CB_FILE_LIST feDataModule feDataRoutines )

foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )
foreach ( FILE ${CB_FILE_LIST} )
  list ( APPEND CB_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${CB_FILE_LIST} )

add_library ( ${LIB_ID} ${F90_SOURCE_FILES} )
if ( GSF_LIBRARY )
  set ( GSF_CB GSF_FEinterface ) # Pre-defined name of call-back library
  add_library ( ${GSF_CB} SHARED ${CB_SOURCE_FILES} )
  target_link_libraries ( ${GSF_CB} vpmCommon_F90 )
  target_link_libraries ( ${LIB_ID} ${GSF_LIBRARY} ${GSF_CB} )
else ( GSF_LIBRARY )
# The module files of the dummy library should not be exported.
# Search for them in the build directory of this library instead.
  include_directories ( BEFORE ${PROJECT_BINARY_DIR}/src/vpmCommon/vpmLinAlg/vpmLinAlgDummy )
  target_link_libraries ( ${LIB_ID} vpmLinAlg_dummy )
endif ( GSF_LIBRARY )
target_link_libraries ( ${LIB_ID} vpmCommon_F90 ${LAPACK_LIBRARIES} )
if ( TARGET SPR_I8 )
  target_link_libraries ( ${LIB_ID} SPR_I8 )
else ( TARGET SPR_I8 )
  target_link_libraries ( ${LIB_ID} SPR )
endif ( TARGET SPR_I8 )
