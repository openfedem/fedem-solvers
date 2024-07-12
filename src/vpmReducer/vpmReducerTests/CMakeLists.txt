# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

if ( BUILD_SOLVER_AS_DLL )
  set ( REDUCER_LIB ${DLL_ID} )
else ( BUILD_SOLVER_AS_DLL )
  set ( REDUCER_LIB ${LIB_ID} )
endif ( BUILD_SOLVER_AS_DLL )
set ( LIB_ID vpmReducerTests )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( USE_MEMPOOL )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_MEMPOOL" )
endif ( USE_MEMPOOL )

include_directories ( ${CMAKE_CURRENT_BINARY_DIR}/.. )

add_library ( ${LIB_ID} smallReducer.f90 solveSupel.f90 )
add_executable ( test_reducer testReducer.C testModels.C ../releaseSingeltons.C )
add_executable ( test_shells testShells.C testModels.C )

target_link_libraries ( ${LIB_ID} vpmReducer_F )
target_link_libraries ( test_shells vpmReducer_F )
target_link_libraries ( test_reducer ${LIB_ID} vpmCommon )
add_dependencies ( check test_reducer )

# Add some reducer regression tests
add_test ( Reduce-Cantilever test_reducer -gvec 0.0 9.81 -dispfile fedem_V.fmx )
add_test ( SmallR-Cantilever test_reducer -small -gvec 0.0 9.81 -resfile smallr.res )
add_test ( Solver-Cantilever test_reducer -solve -gvec 0.0 9.81 -resfile solver.res )
add_test ( BeamPin1 test_reducer 2 -P2 -1 -resfile beampin1.res )
add_test ( BeamPin2 test_reducer 2 -P2 2 -resfile beampin2.res )
add_test ( Rigid-Cylinder test_reducer -P3 2.0 1.0 0.01 2 8 -gvec 0.0 9.81 -resfile cylinder-rgd.res )
add_test ( Flexible-Cylinder test_reducer -P4 2.0 1.0 0.01 2 8 -gvec 0.0 9.81 -resfile cylinder-flx.res )
add_test ( Reduce-Shell1 test_reducer -P5 5.0 1.0 0.01 6 2 -gvec 0.0 0.0 9.81
           -dispfile fedem_V.fmx -resfile shl1-compare.res )
add_test ( Reduce-Shell2 test_reducer -P5 5.0 1.0 0.01 3 1 -gvec 0.0 0.0 9.81 -useParabolic
           -dispfile fedem_V.fmx -resfile shl2-compare.res )
add_test ( Reduce-QuartCyl1 test_reducer -P7 2.0 1.0 0.01 4 5 -gvec 0.0 0.0 9.81
           -dispfile fedem_V.fmx -resfile quartcyl1-compare.res )
add_test ( Reduce-QuartCyl2 test_reducer -P7 2.0 1.0 0.01 2 3 -gvec 0.0 0.0 9.81 -useParabolic
           -dispfile fedem_V.fmx -resfile quartcyl2-compare.res )

file ( GLOB NASTRAN_FILES FEparts/C*_.nas )
foreach ( FILE_NAME ${NASTRAN_FILES} )
  get_filename_component ( MODEL_NAME ${FILE_NAME} NAME_WE )
  configure_file ( ${FILE_NAME} ${MODEL_NAME}.nas COPYONLY )
  file ( STRINGS FEparts/${MODEL_NAME}.sol SOLUTION )
  add_test ( Reduce-${MODEL_NAME} test_reducer -linkfile ${MODEL_NAME}.nas
             -gvec 0.0 9.81 -solution ${SOLUTION} )
endforeach ( FILE_NAME ${NASTRAN_FILES} )

if ( GTest_FOUND )
  # Add some unit tests using gtest (executed via ctest)
  add_executable ( test_FEsolve test_FEsolve.C )
  add_cpp_test ( test_FEsolve ${REDUCER_LIB} )
endif ( GTest_FOUND )

if ( pFUnit_FOUND )
  # Add some unit tests using pFUnit (executed via ctest)
  enable_fortran_tests ()
  add_fortran_test ( test_scatgath testScatGath.pf SAM )
endif ( pFUnit_FOUND )
