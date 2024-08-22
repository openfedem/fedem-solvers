# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID FemLib )
set ( LIB_ID_LIST )
if ( BUILD_TESTS )
  set ( LIB_ID_LIST FemLibTests )
endif ( BUILD_TESTS )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( WIN OR USE_INTEL_FORTRAN )
  string ( APPEND CMAKE_Fortran_FLAGS " -Qdiag-disable:8290,8291" )
else ( WIN OR USE_INTEL_FORTRAN )
  string ( APPEND CMAKE_Fortran_FLAGS " -std=legacy" )
endif ( WIN OR USE_INTEL_FORTRAN )

foreach ( FOLDER ${LIB_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${LIB_ID_LIST} )

if ( USE_INTEL_FORTRAN )
  set ( BLA_VENDOR Intel10_64lp )
endif ( USE_INTEL_FORTRAN )
find_package ( LAPACK REQUIRED )


set ( F90_FILE_LIST pmatStiff cstetra congruenceTrans56 ipri6 isoMatModule
                    andes3ShellModule andes4ShellModule compositeTypeModule )

set ( F77_FILE_LIST beamaux beam dadcos dinv12 fqs ftsa fts hexa hlst
                    ihex ipri itet jaci31 nyteba qbek qmrf tmrf scqs scts )


# Handling f90 compilation
foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )

# Handling f77 compilation
foreach ( FILE ${F77_FILE_LIST} )
  list ( APPEND F77_SOURCE_FILES ${FILE}.f )
endforeach ( FILE ${F77_FILE_LIST} )

add_library ( ${LIB_ID} ${F90_SOURCE_FILES} ${F77_SOURCE_FILES} )

set ( DEPENDENCY_LIST vpmUtilities_F90 SAM ${LAPACK_LIBRARIES} )
target_link_libraries ( ${LIB_ID} ${DEPENDENCY_LIST} )
