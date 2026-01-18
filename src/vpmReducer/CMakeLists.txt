# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( APPL fedem_reducer )
set ( APP2 fedem_partsol )
set ( LIB_ID vpmReducer )
set ( DLL_ID fedem_reducer_core )
set ( LIB_ID_LIST )
if ( BUILD_TESTS )
  set ( LIB_ID_LIST vpmReducerTests )
endif ( BUILD_TESTS )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

option ( BUILD_FMXWRITER_DLL "Build fmx-writer shared library" OFF )
if ( BUILD_FMXWRITER_DLL )
  add_library ( fmxWriter SHARED fmxWriter.f90 )
  target_link_libraries ( fmxWriter LINK_PRIVATE vpmUtilities_F90 )
endif ( BUILD_FMXWRITER_DLL )

if ( BUILD_SOLVER_AS_DLL )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_DLL" )
endif ( BUILD_SOLVER_AS_DLL )

foreach ( FOLDER ${LIB_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${LIB_ID_LIST} )


## Files with header and source with same name
set ( COMPONENT_FILE_LIST )
if ( BUILD_SOLVER_AS_DLL )
  set ( COMPONENT_FILE_LIST ../vpmCommon/FFCmdLineArgInterface )
endif ( BUILD_SOLVER_AS_DLL )

## Pure header files, i.e., header files without a corresponding source file
set ( HEADER_FILE_LIST )

## Pure source files, i.e., without a corresponding header file
set ( SOURCE_FILE_LIST reducerInterface releaseSingeltons )

## Fortran 90 source files
set ( F90_FILE_LIST  cfemFedemModule
                     cfemReducer
                     inputReducerModule
                     rigidModule
                     samReducerModule
                     saveReducerModule
                     wavgmModule
                     reducer
                     solver
                     cmstrs
                     inadd
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

foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )


# Build and install

set ( DEPENDENCY_LIST vpmLinAlg vpmCommon_F90 vpmUtilities_F90 FemLib
                      FFlLib_F FFaLib_F FFaCmdLineArg FFaAlgebra )

add_library ( ${LIB_ID}_F ${F90_SOURCE_FILES} )
target_link_libraries ( ${LIB_ID}_F ${DEPENDENCY_LIST} )
if ( BUILD_SOLVER_AS_DLL )
  add_library ( ${DLL_ID} SHARED ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
  target_link_libraries ( ${DLL_ID} LINK_PRIVATE ${LIB_ID}_F vpmCommon )
else ( BUILD_SOLVER_AS_DLL )
  add_library ( ${LIB_ID} STATIC ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
  target_link_libraries ( ${LIB_ID} LINK_PUBLIC ${LIB_ID}_F vpmCommon )
endif ( BUILD_SOLVER_AS_DLL )

add_executable ( ${APPL} reducermain.C )
add_executable ( ${APP2} solvermain.C )
add_dependencies ( all_solvers ${APPL} ${APP2} )
set ( TARGET_LIST ${APPL} ${APP2} )
if ( BUILD_SOLVER_AS_DLL )
  list ( APPEND TARGET_LIST ${DLL_ID} )
  target_link_libraries ( ${APPL} ${DLL_ID} )
  target_link_libraries ( ${APP2} ${DLL_ID} )
else ( BUILD_SOLVER_AS_DLL )
  target_link_libraries ( ${APPL} ${LIB_ID} )
  target_link_libraries ( ${APP2} ${LIB_ID} )
endif ( BUILD_SOLVER_AS_DLL )
if ( BUILD_FMXWRITER_DLL )
  list ( APPEND TARGET_LIST fmxWriter )
  add_dependencies ( all_solvers fmxWriter )
endif ( BUILD_FMXWRITER_DLL )

install ( TARGETS ${TARGET_LIST}
          RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
          LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" )
