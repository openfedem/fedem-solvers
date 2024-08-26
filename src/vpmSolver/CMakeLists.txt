# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( APPL fedem_solver )
set ( LIB_ID vpmSolver )
set ( DLL_ID fedem_solver_core )
set ( LIB_ID_LIST Control ResponsePos )
if ( BUILD_TEST_REPORTS OR BUILD_TESTS )
  list ( APPEND LIB_ID_LIST vpmSolverTests )
endif ( BUILD_TEST_REPORTS OR BUILD_TESTS )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

option ( USE_AERODYN_PLUGIN "Use the AeroDyn wind model library (NREL)" OFF )
option ( USE_NEMOH_PLUGIN "Use the Nemoh hydrodynamics library (EC Nantes)" OFF )
option ( USE_TNO_PLUGIN "Use the MF-Tyre/Swift tire models (TNO)" OFF )
option ( USE_CTI_PLUGIN "Use the FTire tire model (COSIN)" OFF )
option ( USE_JDTYRE "Use the John Deere tire model" OFF )
option ( USE_PYPLOT "Use the python plot module" OFF )
option ( USE_CONCURRENT_RECOVERY "Build with stress recovery in the loop" OFF )
option ( USE_SP_RECOVERY "Do stress/gage recovery in single precision" OFF )
option ( BUILD_DUMMY_SUPEL_DLL "Build dummy superelement shared library" OFF )
option ( BUILD_CONTROL_AS_DLL "Build the control system module as a DLL" ON )
mark_as_advanced ( USE_AERODYN_PLUGIN USE_NEMOH_PLUGIN USE_TNO_PLUGIN USE_CTI_PLUGIN USE_JDTYRE )
mark_as_advanced ( USE_PYPLOT USE_SP_RECOVERY BUILD_DUMMY_SUPEL_DLL BUILD_CONTROL_AS_DLL )

#
# Find 3rd party packages
#
if ( USE_INTEL_FORTRAN )
  set ( BLA_VENDOR Intel10_64lp )
endif ( USE_INTEL_FORTRAN )
find_package ( BLAS )

if ( USE_AERODYN_PLUGIN )
  find_package ( AeroDyn )
else ( USE_AERODYN_PLUGIN )
  unset ( AeroDyn_LIBRARY CACHE )
  unset ( AeroDyn_DLL CACHE )
  unset ( TURBSIM CACHE )
endif ( USE_AERODYN_PLUGIN )

if ( USE_NEMOH_PLUGIN )
  find_package ( Nemoh )
else ( USE_NEMOH_PLUGIN )
  unset ( Nemoh_LIBRARY CACHE )
endif ( USE_NEMOH_PLUGIN )

if ( USE_TNO_PLUGIN )
  find_package ( TNOTyre )
else ( USE_TNO_PLUGIN )
  unset ( TNO_LIBRARY CACHE )
endif ( USE_TNO_PLUGIN )

if ( USE_CTI_PLUGIN )
  find_package ( FTire )
else ( USE_CTI_PLUGIN )
  unset ( COSIN_LIBRARY CACHE )
endif ( USE_CTI_PLUGIN )

if ( USE_FFTPACK AND TARGET fftpack51 )
  set ( FFT_LIBRARY fftpack51 )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_FFTPACK" )
else ( USE_FFTPACK AND TARGET fftpack51 )
  message ( STATUS "Configuring without frequency domain analysis support" )
endif( USE_FFTPACK AND TARGET fftpack51 )

if ( WIN )
  string ( APPEND CMAKE_EXE_LINKER_FLAGS " /NODEFAULTLIB:LIBCMT" )
endif ( WIN )
if ( Nemoh_LIBRARY AND WIN )
  string ( APPEND CMAKE_EXE_LINKER_FLAGS " /STACK:4000000" )
endif ( Nemoh_LIBRARY AND WIN )
if ( BLAS_LIBRARIES MATCHES mkl_intel )
  message ( STATUS "Using Fast Fourier Transformation implementation from the Intel Math Kernal Library" )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_MKL" )
endif ( BLAS_LIBRARIES MATCHES mkl_intel )
if ( USE_JDTYRE )
  list ( APPEND LIB_ID_LIST JDTire )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_JDTYRE" )
endif ( USE_JDTYRE )
if ( USE_CONCURRENT_RECOVERY )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_HAS_RECOVERY" )
  if ( USE_SP_RECOVERY )
    string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_RECOVERY=1" )
  else ( USE_SP_RECOVERY )
    string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_RECOVERY=2" )
  endif ( USE_SP_RECOVERY )
endif ( USE_CONCURRENT_RECOVERY )

foreach ( FOLDER ${LIB_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${LIB_ID_LIST} )


## Files with header and source with same name
set ( COMPONENT_FILE_LIST )
## Pure header files, i.e., header files without a corresponding source file
set ( HEADER_FILE_LIST )
## Pure source files, i.e., without a corresponding header file
set ( SOURCE_FILE_LIST HWAFLSmanager )

## Fortran 90 source files
set ( F90_FILE_LIST  addInSysModule
                     bushingElementRoutinesModule bushingElementTypeModule
                     contactElementRoutinesModule contactElementTypeModule
                     contactSurfaceModule controlRoutinesModule controlTypeModule
                     corotUtilModule curvePointTypeModule curveTypeModule
                     damperRoutinesModule damperTypeModule
                     dbgUnitsModule diffractionModule
                     energyInt engineRoutinesModule environmentTypeModule
                     finiteElementModule forceRoutinesModule forceTypeModule
                     freqResponseModule frictionRoutinesModule frictionTypeModule
                     functionTypeModule genericPartModule hydroDynamicsModule
                     initiateBushingElmTypeModule initiateContactElmTypeModule
                     initiateDamperTypeModule initiateFrictionTypeModule
                     initiateFunctionTypeModule initiateJointTypeModule
                     initiateModule initiateModesTypeModule initiateRoadTypeModule
                     initiateSensorTypeModule initiateSpringTypeModule
                     initiateSupElTypeModule initiateSystemTypeModule
                     initiateTireTypeModule initiateTriadTypeModule
                     initiateUserdefElTypeModule initiateWindTurbineModule
                     inverseModule massMatCorrModule
                     massRoutinesModule massTypeModule
                     masterSlaveJointRoutinesModule masterSlaveJointTypeModule
                     mechanismTypeModule modesRoutinesModule modesTypeModule
                     motionRoutinesModule motionTypeModule
                     newmarkRoutinesModule
                     normRoutinesModule normTypeModule
                     profilerModule restartModule
                     roadRoutinesModule roadTypeModule
                     samSolverModule saveModule saveNastranModule saveVTFModule
                     sensorTypeModule
                     solverDriver solverModule solverRoutinesModule
                     springRoutinesModule springTypeModule
                     staticEquilibriumModule
                     supElLoadTypeModule supElRoutinesModule
                     systemTypeModule timeStepModule
                     tireRoutinesModule tireTypeModule
                     userdefElRoutinesModule userdefElTypeModule
                     windTurbineRoutinesModule windTurbineTypeModule
                     FNVmodule FNVwaveForceModule
    )

if ( USE_EXT_CTRLSYS )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_EXTCTRL" )
  list ( APPEND F90_FILE_LIST extCtrlSysRoutinesModule extCtrlSysTypeModule )
endif ( USE_EXT_CTRLSYS )

if ( USE_JDTYRE )
  list ( APPEND F90_FILE_LIST JDTire/JD_tireRoutines JDTire/JD_tireRoutinesModule )
else ( USE_JDTYRE )
  list ( APPEND F90_FILE_LIST JDTire/JD_dummy )
endif ( USE_JDTYRE )

if ( USE_PYPLOT )
  list ( APPEND F90_FILE_LIST pyplot_module )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_PYPLOT" )
else ( USE_PYPLOT )
  list ( APPEND F90_FILE_LIST dummy_pyplot )
endif ( USE_PYPLOT )

if ( NOT ( BUILD_SOLVER_AS_DLL AND BUILD_DUMMY_SUPEL_DLL ) )
  list ( APPEND F90_FILE_LIST dummyLoadSupelMatrix )
endif ( NOT ( BUILD_SOLVER_AS_DLL AND BUILD_DUMMY_SUPEL_DLL ) )

if ( USE_CONCURRENT_RECOVERY )
  set ( REC_FILE_LIST stressRecoveryModule
                     ../vpmStress/displacementModule
                     ../vpmStress/samStressModule
                     ../vpmStress/strainAndStressUtils
                     ../vpmStress/stressRoutines
                     ../vpmStress/elStressModule
                     ../vpmStress/saveStressModule
                     ../vpmStress/saveStrainGageModule
                     ../vpmStress/saveVTFModule
                     ../vpmStress/strainRosetteModule
                     dummyResStressModule
    )
endif ( USE_CONCURRENT_RECOVERY )


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
foreach ( FILE ${REC_FILE_LIST} )
  list ( APPEND REC_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )


if ( AeroDyn_LIBRARY )
  list ( APPEND EXTERNAL_LIBRARIES ${AeroDyn_LIBRARY} )
endif ( AeroDyn_LIBRARY )
if ( Nemoh_LIBRARY )
  list ( APPEND EXTERNAL_LIBRARIES ${Nemoh_LIBRARY} )
endif ( Nemoh_LIBRARY )
if ( COSIN_LIBRARY )
  list ( APPEND EXTERNAL_LIBRARIES ${COSIN_LIBRARY} )
endif ( COSIN_LIBRARY )
if ( TNO_LIBRARY )
  list ( APPEND EXTERNAL_LIBRARIES ${TNO_LIBRARY} )
endif ( TNO_LIBRARY )
if ( VTF_LIBRARY )
  list ( APPEND EXTERNAL_LIBRARIES ${VTF_LIBRARY} )
endif ( VTF_LIBRARY )
if ( FFT_LIBRARY )
  list ( APPEND EXTERNAL_LIBRARIES ${FFT_LIBRARY} )
endif ( FFT_LIBRARY )

set ( FT_COMMON_LIBRARIES FFrLib_F )
set ( FT_KERNEL_LIBRARIES vpmLinAlg
                          FFaFunctionLib_F90
                          FiUserElmPlugin_F90
    )
set ( REC_COMMONLIBS FemLib FFlLib_F )
set ( REC_KERNELLIBS vpmCommon_F90 FiDeviceFunctions_F90 FFpLib_F )

if ( USE_CONCURRENT_RECOVERY )
  if ( LINUX )
    # Add stress recovery routines as a separate library
    # to resolve dependency problem
    add_library ( ${LIB_ID}_Rec ${REC_SOURCE_FILES} )
    target_link_libraries ( ${LIB_ID}_Rec ${REC_KERNELLIBS} ${REC_COMMONLIBS} )
    list ( INSERT FT_KERNEL_LIBRARIES 0 ${LIB_ID}_Rec )
  else ( LINUX )
    list ( APPEND F90_SOURCE_FILES ${REC_SOURCE_FILES} )
    list ( INSERT FT_COMMON_LIBRARIES 0 ${REC_COMMONLIBS} )
    list ( INSERT FT_KERNEL_LIBRARIES 0 ${REC_KERNELLIBS} )
  endif ( LINUX )
else ( USE_CONCURRENT_RECOVERY )
  list ( INSERT FT_COMMON_LIBRARIES 0 FemLib )
  list ( INSERT FT_KERNEL_LIBRARIES 0 ${REC_KERNELLIBS} )
endif ( USE_CONCURRENT_RECOVERY )

if ( BUILD_CONTROL_AS_DLL )
  set ( FT_LIBRARIES fedem_solver_ctrl )
else ( BUILD_CONTROL_AS_DLL )
  set ( FT_LIBRARIES Control )
endif ( BUILD_CONTROL_AS_DLL )
list ( APPEND FT_LIBRARIES ${FT_KERNEL_LIBRARIES} ${FT_COMMON_LIBRARIES} )
set ( DEPENDENCY_LIST ${FT_LIBRARIES} ${EXTERNAL_LIBRARIES} )
if ( USE_JDTYRE )
  list ( INSERT DEPENDENCY_LIST 0 JDTire_F )
endif ( USE_JDTYRE )


# Build and install

add_library ( ${LIB_ID}_F ${F90_SOURCE_FILES} )
add_library ( ${LIB_ID}_C ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
if ( BUILD_SOLVER_AS_DLL )
  add_library ( ${DLL_ID} SHARED solverInterface.C solverInterface.h )
  if ( BUILD_DUMMY_SUPEL_DLL )
    add_library ( vpmSupEls SHARED dummyLoadSupelMatrix.f90 )
    list ( INSERT DEPENDENCY_LIST 0 vpmSupEls )
  endif ( BUILD_DUMMY_SUPEL_DLL )
else ( BUILD_SOLVER_AS_DLL )
  add_library ( ${LIB_ID} STATIC solverInterface.C solverInterface.h )
endif ( BUILD_SOLVER_AS_DLL )

target_link_libraries ( ${LIB_ID}_F ${DEPENDENCY_LIST} )
if ( BUILD_SOLVER_AS_DLL )
  target_link_libraries ( ${DLL_ID} LINK_PRIVATE ${LIB_ID}_F ${LIB_ID}_C vpmCommon )
else ( BUILD_SOLVER_AS_DLL )
  target_link_libraries ( ${LIB_ID} LINK_PUBLIC ${LIB_ID}_F ${LIB_ID}_C vpmCommon )
endif ( BUILD_SOLVER_AS_DLL )

add_executable ( ${APPL} main.C vpmSolverTests/compareResponse.C )
add_dependencies ( all_solvers ${APPL} )
set ( TARGET_LIST ${APPL} )
if ( BUILD_SOLVER_AS_DLL )
  target_link_libraries ( ${APPL} ${DLL_ID} )
  list ( APPEND TARGET_LIST ${DLL_ID} )
  if ( BUILD_DUMMY_SUPEL_DLL )
    list ( APPEND TARGET_LIST vpmSupEls )
  endif ( BUILD_DUMMY_SUPEL_DLL )
else ( BUILD_SOLVER_AS_DLL )
  target_link_libraries ( ${APPL} ${LIB_ID} )
endif ( BUILD_SOLVER_AS_DLL )

install ( TARGETS ${TARGET_LIST}
          RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
          LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" )

if ( AeroDyn_DLL )
  install ( FILES ${AeroDyn_DLL}
            DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" )
elseif ( AeroDyn_LIBRARY )
  install ( FILES ${AeroDyn_LIBRARY}
            DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" )
endif ( AeroDyn_DLL )
if ( TURBSIM )
  install ( PROGRAMS ${TURBSIM}
            DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" )
  get_filename_component ( TURBSIM_DIR "${TURBSIM}" PATH )
  get_filename_component ( TURBSIM_DIR "${TURBSIM_DIR}" PATH )
  file ( GLOB TEMPLATES ${TURBSIM_DIR}/Templates/*.inp )
  install ( FILES ${TEMPLATES}
            DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/Templates" )
endif ( TURBSIM )
