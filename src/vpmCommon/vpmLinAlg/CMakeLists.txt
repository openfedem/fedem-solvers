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

if ( USE_LINALG_GSF )
  find_package ( FTLinAlg )
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
if ( NOT FTLINALG_FOUND )
# When the 3rd-party equation solver is not available,
# add a dummy library containing the expected symbols.
# This library is compiled in a sub-folder, to avoid that
# its generated module files are mixed with the others.
  list ( INSERT SUBFOLDER_LIST 0 vpmLinAlgDummy )
endif ( NOT FTLINALG_FOUND )
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
if ( FTLINALG_FOUND )
  add_library ( ${LIB_ID}_CB ${CB_SOURCE_FILES} )
  target_link_libraries ( ${LIB_ID}_CB vpmCommon_F90 )
  target_link_libraries ( ${LIB_ID} ${FTLINALG_LIBRARIES} ${LIB_ID}_CB )
else ( FTLINALG_FOUND )
# The module files of the dummy library should not be exported.
# Search for them in the build directory of this library instead.
  include_directories ( BEFORE ${PROJECT_BINARY_DIR}/src/vpmCommon/vpmLinAlg/vpmLinAlgDummy )
  target_link_libraries ( ${LIB_ID} vpmLinAlg_dummy )
endif ( FTLINALG_FOUND )
target_link_libraries ( ${LIB_ID} vpmCommon_F90 ${LAPACK_LIBRARIES} )
if ( TARGET SPR_I8 )
  target_link_libraries ( ${LIB_ID} SPR_I8 )
else ( TARGET SPR_I8 )
  target_link_libraries ( ${LIB_ID} SPR )
endif ( TARGET SPR_I8 )
