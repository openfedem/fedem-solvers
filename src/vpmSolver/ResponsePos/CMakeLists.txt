# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( LIB_ID ResponsePos )
set ( APP_ID response_pos )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )


## Fortran 90 source files
set ( F90_FILE_LIST  response_pos
                     ../bushingElementTypeModule
                     ../contactElementTypeModule
                     ../contactSurfaceModule
                     ../controlTypeModule
                     ../curvePointTypeModule
                     ../damperTypeModule
                     ../environmentTypeModule
                     ../forceTypeModule
                     ../frictionTypeModule
                     ../functionTypeModule
                     ../massTypeModule
                     ../masterSlaveJointTypeModule
                     ../mechanismTypeModule
                     ../modesTypeModule
                     ../motionTypeModule
                     ../normTypeModule
                     ../profilerModule
                     ../roadTypeModule
                     ../saveModule
                     ../sensorTypeModule
                     ../springTypeModule
                     ../supElLoadTypeModule
                     ../systemTypeModule
                     ../tireTypeModule
                     ../userdefElTypeModule
                     ../windTurbineTypeModule
                     ../../vpmStress/saveStrainGageModule
                     ../../vpmStress/strainRosetteModule
                     ../../vpmStress/strainAndStressUtils
                     dummies
    )


foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )


# Build and install

add_library ( ${LIB_ID} ${F90_SOURCE_FILES} )
add_executable ( ${APP_ID} main.C )
target_link_libraries ( ${LIB_ID} vpmCommon_F90 FFaFunctionLib_F90 )
target_link_libraries ( ${APP_ID} ${LIB_ID} )

add_custom_command ( OUTPUT response_pos.frs
                     COMMAND ${APP_ID} DEPENDS ${APP_ID} )
add_custom_target ( response_pos_frs ALL
                    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/response_pos.frs )

install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/response_pos.frs
          DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/resources" )
