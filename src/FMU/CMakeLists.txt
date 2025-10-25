# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID fedem_fmu )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

string ( APPEND CMAKE_CXX_FLAGS_DEBUG " -DFMU_DEBUG" )
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL GNU )
  string ( APPEND CMAKE_CXX_FLAGS " -Wno-unused-parameter" )
endif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL GNU )

# Build and install

include_directories ( "${CMAKE_CURRENT_SOURCE_DIR}/FMI2/headers" )
add_library ( ${LIB_ID} SHARED fmu_template.C )
target_link_libraries ( ${LIB_ID} Admin )

set ( CLOUD_DIR "${CMAKE_INSTALL_PREFIX}/bin/Templates/cloudsim" )
install ( TARGETS ${LIB_ID}
          RUNTIME DESTINATION "${CLOUD_DIR}"
          LIBRARY DESTINATION "${CLOUD_DIR}" )
