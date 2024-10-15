# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

################################################################################
# This is the top-level cmake project file for the Fedem FMU module.
################################################################################

cmake_minimum_required ( VERSION 2.8...3.5 )
if ( POLICY CMP0076 )
  cmake_policy ( SET CMP0076 NEW ) # convert relative target source path names
endif ( POLICY CMP0076 )

# Project setup

set ( APPLICATION_ID fedemFMU )
set ( DOMAIN_ID FEDEM )
set ( PACKAGE_ID SOLVERS )

project ( ${APPLICATION_ID} CXX )
message ( STATUS "Generating build project for ${PROJECT_SOURCE_DIR}" )
get_filename_component ( PARENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY )

find_path ( _MODULES FedemConfig.cmake
            PATHS $ENV{CMAKE_MODULES}
           "${PARENT_SOURCE_DIR}/fedem-foundation/cmake/Modules/"
          )
if ( _MODULES )
  message ( STATUS "NOTE : Using ${_MODULES}" )
  list ( APPEND CMAKE_MODULE_PATH ${_MODULES} )
else ( _MODULES )
  message ( STATUS "ERROR : Missing path to FedemConfig.cmake" )
  message ( FATAL_ERROR "Set environment variable CMAKE_MODULES and try again" )
endif ( _MODULES )
unset ( _MODULES CACHE )

include ( FedemConfig )

add_subdirectory ( ../src/FMU "${CMAKE_CURRENT_BINARY_DIR}/src" )
