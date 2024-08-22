# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

################################################################################
# This is the top-level cmake project file for the Fedem solver sources.
################################################################################

cmake_minimum_required ( VERSION 2.8...3.5 )
if ( POLICY CMP0076 )
  cmake_policy ( SET CMP0076 NEW ) # convert relative target source path names
endif ( POLICY CMP0076 )

# Project setup

set ( APPLICATION_ID fedemKernel )
set ( DOMAIN_ID FEDEM )
set ( PACKAGE_ID SOLVERS )

project ( ${APPLICATION_ID} CXX C Fortran )
message ( STATUS "Generating build project for ${PROJECT_SOURCE_DIR}" )

find_path ( _MODULES FedemConfig.cmake
            PATHS $ENV{CMAKE_MODULES}
           "${PROJECT_SOURCE_DIR}/fedem-foundation/cmake/Modules/"
          )
if ( _MODULES )
  message ( STATUS "NOTE : Using ${_MODULES}" )
  list ( APPEND CMAKE_MODULE_PATH ${_MODULES} )
else ( _MODULES )
  message ( STATUS "ERROR : Missing path to FedemConfig.cmake" )
  message ( FATAL_ERROR "Set environment variable CMAKE_MODULES and try again" )
endif ( _MODULES )
unset ( _MODULES CACHE )

include ( CheckPFUnit )
include ( FedemConfig )

option ( USE_MEMPOOL "Use memory pool for heap allocation in FE library" OFF )
option ( USE_PROFILER "Use CPU and Memory profiler" OFF )
option ( BUILD_TEST_REPORTS "Enable test and profiling reports" OFF )
mark_as_advanced ( USE_MEMPOOL USE_PROFILER BUILD_TEST_REPORTS )


# Enable unit and regression testing

set ( CTEST_OPTIONS --force-new-ctest-process --output-on-failure -O CTest.txt )
if ( BUILD_TEST_REPORTS )
  list ( APPEND CTEST_OPTIONS --test-action Test --no-compress-output )
  # Needed by the profiling execution (see below)
  set ( MEMORYCHECK_COMMAND_OPTIONS "--tool=callgrind -v" )
  include ( CTest )
elseif ( BUILD_TESTS )
  enable_testing ()
endif ( BUILD_TEST_REPORTS )
if ( CMAKE_CONFIGURATION_TYPES )
  list ( APPEND CTEST_OPTIONS --build-config \"$<CONFIGURATION>\" )
endif ( CMAKE_CONFIGURATION_TYPES )

if ( BUILD_TEST_REPORTS OR BUILD_TESTS )

# The target check is to be used (instead of test) to ensure that
# all test executables are built before executing the tests.
  add_custom_target ( check COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS} )

  include ( GTest ) # Using the google test framework for C++ unit tests
  include ( pFUnit ) # Using the pFUnit test framework for Fortran unit tests

endif ( BUILD_TEST_REPORTS OR BUILD_TESTS )

# Setting up profiling with callgrind.
# Set PROF_RANGE for running a specific number of tests by numbers
# (must be specified at cmake level, see CTest documentation).
# cmake ....... -DPROF_RANGE=3,3 (for example)

if ( BUILD_TEST_REPORTS )
  set ( CTEST_OPTIONS --force-new-ctest-process
                      -D ExperimentalMemCheck -O CTest.log )
  if ( PROF_RANGE )
    list ( APPEND CTEST_OPTIONS -I ${PROF_RANGE} )
  endif ( PROF_RANGE )
  add_custom_target ( code_profile
                      COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS}
                      COMMENT "Generating profiling with callgrind"
                      VERBATIM )
endif ( BUILD_TEST_REPORTS )


# Generate documentation for code coverage

if ( LINUX AND BUILD_TEST_REPORTS )
  include ( CodeCoverage )

# skip coverage reporting for all sources under these folders
  set ( EXCLUDES src/vpmReducer src/vpmStress
                 src/vpmCommon/vpmLinAlg/vpmLinAlgDummy
                 fedem-foundation/src/FemLib fedem-foundation/src/FFlLib
                 fedem-foundation/src/Admin *vpmLinAlg/fe *Tests )

# build html coverage report with gcovr
  gcovr_target ( NAME gcoverage HTML EXCLUDE ${EXCLUDES} )

# build xml (Cobertura) coverage report with gcovr
  gcovr_target ( NAME gcoverage XML EXCLUDE ${EXCLUDES} )

# build xml (Sonarqube) coverage report with gcovr
  gcovr_target ( NAME gcoverage SONAR EXCLUDE ${EXCLUDES} )

# build html and xml (both Coberture and Sonarqube) coverage report with gcovr
  gcovr_target ( NAME gcoverage ALL EXCLUDE ${EXCLUDES} )

# build html report with lcov
  list ( APPEND EXCLUDES src/OP2converter
                         src/vpmGraphExp src/vpmSolver/ResponsePos
                         ${CMAKE_CURRENT_BINARY_DIR}/Admin )
  lcov_html_target ( NAME lcoverage EXCLUDE ${EXCLUDES} /usr/include
                     EXECUTABLE ${CMAKE_CTEST_COMMAND} -j ${PROCESSOR_COUNT} )
endif ( LINUX AND BUILD_TEST_REPORTS )


# Generate source code documentation with doxygen, if available

find_program ( DOXYGEN "doxygen" )
if ( DOXYGEN )
  message ( STATUS "Found Doxygen : ${DOXYGEN}" )
  if ( WIN )
    execute_process ( COMMAND "echo %DATE:~7,2%.%DATE:~4,2%.%DATE:~10,4%"
                      OUTPUT_VARIABLE BUILD_DATE )
  else ( WIN )
    execute_process ( COMMAND "date" "+%e %b %Y"
                      OUTPUT_VARIABLE BUILD_DATE )
  endif ( WIN )
  configure_file ( src/Doxyfile.solver.in Doxyfile.solver )
  configure_file ( src/vpmSolver/mainpage.dox.in mainpage.dox )
  add_custom_target ( doc )
  add_custom_target ( solver_doc
                      COMMAND ${DOXYGEN} ${PROJECT_BINARY_DIR}/Doxyfile.solver
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      COMMENT "Generating solver source code documentation"
                      VERBATIM )
  add_dependencies ( doc solver_doc )
endif ( DOXYGEN )


# Process the sub-folders

set ( USE_FORTRAN true )
set ( SRC_DIR "${CMAKE_CURRENT_BINARY_DIR}/src" )

add_subdirectory ( fedem-foundation/src "${CMAKE_CURRENT_BINARY_DIR}/fedem-foundation" )
add_subdirectory ( SAM/src "${CMAKE_CURRENT_BINARY_DIR}/src/SAM" )
add_subdirectory ( src )
if ( BUILD_TEST_REPORTS OR BUILD_TESTS )
  add_subdirectory ( PythonAPI )
endif ( BUILD_TEST_REPORTS OR BUILD_TESTS )

# Install the Fortran redistributables

if ( WIN AND DEFINED ENV{IFORT_COMPILER23} )
  file ( TO_CMAKE_PATH $ENV{IFORT_COMPILER23}redist/intel64_win/compiler IFORT_PATH )
  file ( GLOB FTN_REDIST ${IFORT_PATH}/*md.dll )
  if ( FTENV_VERBOSE )
    message ( STATUS "Installing Fortran redistributables from ${IFORT_PATH} to ${CMAKE_INSTALL_PREFIX}" )
  endif ( FTENV_VERBOSE )
  install ( FILES ${FTN_REDIST} DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" )
endif ( WIN AND DEFINED ENV{IFORT_COMPILER23} )
if ( WIN AND DEFINED ENV{MKLROOT} )
  file ( TO_CMAKE_PATH $ENV{MKLROOT}/redist/intel64 MKL_PATH )
  if ( FTENV_VERBOSE )
    message ( STATUS "Installing MKL redistributables from ${MKL_PATH} to ${CMAKE_INSTALL_PREFIX}" )
  endif ( FTENV_VERBOSE )
  set ( MKL_REDIST ${MKL_PATH}/mkl_core.2.dll ${MKL_PATH}/mkl_def.2.dll
                   ${MKL_PATH}/mkl_avx2.2.dll ${MKL_PATH}/mkl_vml_avx2.2.dll
                   ${MKL_PATH}/mkl_avx512.2.dll ${MKL_PATH}/mkl_vml_avx512.2.dll
                   ${MKL_PATH}/mkl_intel_thread.2.dll )
  install ( FILES ${MKL_REDIST} DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" )
endif ( WIN AND DEFINED ENV{MKLROOT} )
