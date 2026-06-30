# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

if ( BUILD_SOLVER_AS_DLL )
  set ( SOLVER_LIB ${DLL_ID} )
else ( BUILD_SOLVER_AS_DLL )
  set ( SOLVER_LIB ${LIB_ID} )
endif ( BUILD_SOLVER_AS_DLL )

set ( LIB_ID vpmSolverTests )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

#
# Specialized executables for regression testing (see below)
#
set ( TEST1 test_external )
set ( TEST2 test_restart )
set ( TEST3 test_slviter )
set ( TEST4 test_inverse )
set ( TEST5 test_memchk )
add_executable ( ${TEST1} testExternal.C ../readFSI.C compareResponse.C )
add_executable ( ${TEST2} testRestart.C compareResponse.C )
add_executable ( ${TEST3} testSolveIter.C )
add_executable ( ${TEST4} testInverse.C compareResponse.C xerbla.f )
add_executable ( ${TEST5} testMemcheck.C )
target_link_libraries ( ${TEST1} ${SOLVER_LIB} )
target_link_libraries ( ${TEST2} ${SOLVER_LIB} )
target_link_libraries ( ${TEST3} ${SOLVER_LIB} )
target_link_libraries ( ${TEST4} ${SOLVER_LIB} )
target_link_libraries ( ${TEST5} ${SOLVER_LIB} )

#
# Unit testing
#
if ( pFUnit_FOUND )
# Add some unit tests using pFUnit (executed via ctest)
  enable_fortran_tests ( SRCDIR )
  add_fortran_test ( test_readfsi testFSIreader.pf SRCDIR ${SOLVER_LIB} )
  add_fortran_test ( test_sysbeam testBeam.pf SRCDIR vpmCommon_F90 FemLib )
  add_fortran_test ( test_motion testMotion.pf SRCDIR vpmCommon_F90 )
  target_include_directories ( test_readfsi PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/.. )
  target_sources ( test_sysbeam PUBLIC ../finiteElementModule.f90 dummyModules.f90 )
  target_sources ( test_motion PUBLIC ../motionTypeModule.f90 ../motionRoutinesModule.f90 dummyModules.f90 )
  add_subdirectory ( moreTests )
endif ( pFUnit_FOUND )

if ( GTest_FOUND )
# Add some unit tests using gtest (executed via ctest)
  add_executable ( test_solver test_solver.C )
  add_cpp_test ( test_solver ${SOLVER_LIB} )
endif ( GTest_FOUND )

#
# Regression testing.
# The test cases are maintained in the solverTests repository,
# which is included as a submodule in the ftKernel repository.
#
if ( EXISTS ${PROJECT_SOURCE_DIR}/solverTests/CMakeLists.txt )
  add_subdirectory ( ${PROJECT_SOURCE_DIR}/solverTests
                    "${CMAKE_CURRENT_BINARY_DIR}/solverTests" )

  # The target 'check' ensures that the test programs are up to date, and
  # copies all input files to the build directory before executing the tests
  add_dependencies ( check fedem_reducer fedem_partsol fedem_solver
                     test_external test_restart test_inverse )
  if ( TARGET copy_test_input_files )
    add_dependencies ( check copy_test_input_files )
  endif ( TARGET copy_test_input_files )

endif ( EXISTS ${PROJECT_SOURCE_DIR}/solverTests/CMakeLists.txt )
