# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

set ( PACKAGE_ID PythonAPI )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( CMAKE_VERSION VERSION_LESS 3.12 )
  find_package ( PythonInterp 3 ) # sets ${PYTHON_EXECUTABLE}
else ( CMAKE_VERSION VERSION_LESS 3.12 )
  find_package ( Python3 ) # sets ${Python3_EXECUTABLE}
  set ( PythonInterp_FOUND ${Python3_FOUND} )
  set ( PYTHON_EXECUTABLE ${Python3_EXECUTABLE} )
endif ( CMAKE_VERSION VERSION_LESS 3.12 )

if ( PythonInterp_FOUND )
  set ( PYTHON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )
  add_subdirectory ( PythonAPITests )
endif ( PythonInterp_FOUND )
