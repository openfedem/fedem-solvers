# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID vpmUtilities )
set ( LIB_ID_LIST )
if ( BUILD_TESTS )
  set ( LIB_ID_LIST vpmUtilitiesTests )
endif ( BUILD_TESTS )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

foreach ( FOLDER ${LIB_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${LIB_ID_LIST} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

string ( APPEND CMAKE_CXX_FLAGS " -D_MAX_DBFIL=256" )
if ( WIN )
  string ( APPEND CMAKE_Fortran_FLAGS " -Qdiag-disable:8290" )
endif ( WIN )


## Files with header and source with same name
set ( COMPONENT_FILE_LIST )

## Pure header files, i.e., header files without a corresponding source file
set ( HEADER_FILE_LIST )

## Pure implementation files, i.e., source files without corresponding header
set ( SOURCE_FILE_LIST computerConfig sparseMatrixUtils )

## Special C files
set ( H_FILE_LIST )
set ( C_FILE_LIST binaryDB clksec )

## Fortran 90 source files
set ( F90_FILE_LIST kindModule
                    rotationModule
                    binaryDBInterface
                    fileUtilitiesModule
                    manipMatrixModule
                    reportErrorModule
                    computerConfigInterface
                    allocationModule
                    timerModule
                    scratchArrayModule
                    searchAndSortModule
                    sparseMatrixInterface
                    sparseMatrixModule
                    denseMatrixModule
                    progressModule
                    diskMatrixModule
                    sDiskMatrixModule
                    inputUtilities
                    closeLogFiles
                    wrimsg
    )

foreach ( FILE ${COMPONENT_FILE_LIST} )
  list ( APPEND CPP_SOURCE_FILES ${FILE}.C )
  list ( APPEND CPP_HEADER_FILES ${FILE}.H )
endforeach ( FILE ${COMPONENT_FILE_LIST} )

foreach ( FILE ${HEADER_FILE_LIST} )
  list ( APPEND CPP_HEADER_FILES ${FILE}.H )
endforeach ( FILE ${HEADER_FILE_LIST} )

foreach ( FILE ${SOURCE_FILE_LIST} )
  list ( APPEND CPP_SOURCE_FILES ${FILE}.C )
endforeach ( FILE ${SOURCE_FILE_LIST} )

foreach ( FILE ${H_FILE_LIST} )
  list ( APPEND CPP_HEADER_FILES ${FILE}.h )
endforeach ( FILE ${H_FILE_LIST} )

foreach ( FILE ${C_FILE_LIST} )
  list ( APPEND CPP_SOURCE_FILES ${FILE}.c )
endforeach ( FILE ${C_FILE_LIST} )

foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )


add_library ( ${LIB_ID} ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
add_library ( ${LIB_ID}_F90 ${F90_SOURCE_FILES} )
target_link_libraries ( ${LIB_ID} FFaOS )
target_link_libraries ( ${LIB_ID}_F90 ${LIB_ID} FFaLib_F SAM )
