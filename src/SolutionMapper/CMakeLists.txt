# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( APPL fedem_solmap )
set ( LIB_ID SolMapper )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

set ( DEPENDENCY_LIST FFrLib FFlIOAdaptors FFaCmdLineArg FFaAlgebra Admin )

if ( USE_MEMPOOL )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_MEMPOOL" )
endif ( USE_MEMPOOL )

# Need to use FFlLib classes in the FTK namespace
string ( APPEND CMAKE_CXX_FLAGS " -DFF_NAMESPACE=FTK" )

string ( APPEND CMAKE_CXX_FLAGS_DEBUG " -DFT_DEBUG=1" )


# Build and install

add_executable ( ${APPL} main.C ../vpmCommon/cmdLineArgInitStd.C )
target_link_libraries ( ${APPL} ${DEPENDENCY_LIST} )
add_dependencies ( all_solvers ${APPL} )
install ( TARGETS ${APPL} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" )
