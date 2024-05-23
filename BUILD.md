<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# FEDEM Solvers Build

Currently, we support building on Windows
using Microsoft Visual Studio 2019 (or later) and Intel&reg; Fortran Compilers,
and on Linux using the GNU compilers (gcc and gfortran).
The build system is based on the cross-platform [CMake](https://cmake.org/) tool.
On Windows, this tool is embedded in the Visual Studio 2019 installation.
On Linux, you can install `cmake` from the package manager of your distribution,
e.g., on Ubuntu

    sudo apt install cmake

## Build instruction

The top level `CMakeLists.txt` file of this repository contains
the necessary setup and configuration. It also controls the execution of
the unit- and regression test [see below](#regression-testing).

We use the packages [googletest](https://github.com/google/googletest) and
[pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit) to implement
some unit tests for the C++ and Fortran code, respectively. Therefore,
you need to set the environment variables `GTEST_ROOT` and `PFUNIT` to point to
the root path of valid installations of the googletest and pFUnit frameworks,
respectively, before executing the `cmake` commands shown below. See
[Installing googletest and pFUnit](https://github.com/openfedem/fedem-foundation/blob/main/pFUnit/README.md#installing-googletest-and-pfunit)
for instructions on how to install these packages on your system.

To build the FEDEM solvers from the sources, proceed as follows:

- From a bash shell or command prompt, clone the sources of this repository:

      mkdir ~/Fedem-src
      cd ~/Fedem-src
      git clone --recurse-submodules git@github.com:openfedem/fedem-solvers.git

- On Linux, the FEDEM solvers and tests can be built and executed by:

      cd ~/Fedem-src/fedem-solvers
      mkdir Release; cd Release
      cmake .. -DCMAKE_BUILD_TYPE=Release \
               -DCMAKE_INSTALL_PREFIX=$HOME \
               -DBUILD_SOLVER_AS_DLL=ON -DBUILD_CONTROL_AS_DLL=ON -DBUILD_TEST_REPORTS=ON \
               -DUSE_CONCURRENT_RECOVERY=ON -DUSE_SP_RECOVERY=ON -DUSE_FFTPACK=ON \
               -DFT_LARGE_MODELS=ON -DFT_TOLERANCE=1.0e-12
      make check

  Then you can run `make install` to install the solver executables in the `$HOME/bin` folder.
  The shared object libraries (`*.so` files) will be installed under `$HOME/lib`.
  If a different installation location is wanted, you need to change the
  `-DCMAKE_INSTALL_PREFIX=$HOME` line in the `cmake` command above.

- On Windows, use the following bat script to configure the build with
  Visual Studio 2019 and the Intel&reg; Fortran Compilers.

      @echo off
      title FEDEM solvers configuration
      call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019
      set /p VERSION=<%USERPROFILE%\Fedem-src\fedem-solvers\cfg\VERSION
      "%VS2019INSTALLDIR%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" ^
      -G "Visual Studio 16 2019" ^
      -S %USERPROFILE%\Fedem-src\fedem-solvers ^
      -B %USERPROFILE%\Fedem-build\solvers ^
      -DCMAKE_INSTALL_PREFIX=%USERPROFILE%\Fedem-install\%VERSION% ^
      -DBUILD_SOLVER_AS_DLL=ON -DBUILD_CONTROL_AS_DLL=ON ^
      -DUSE_CONCURRENT_RECOVERY=ON -DUSE_SP_RECOVERY=ON -DUSE_FFT_PACK=ON ^
      -DFT_LARGE_MODELS=ON -DFT_TOLERANCE=1.0e-10
      pause

- Then (still on Windows), open the generated solution file
  `%USERPROFILE%\Fedem-build\solvers\fedemKernel.sln`
  in Visual Studio and build the `all_solvers` target for `Release`
  configuration to compile all solver modules. Build the `check` target
  if you also want to execute the tests after building the solvers.

  Build the `INSTALL` target to install the resulting binaries
  (`.exe` file and dependent `.dll` files) in the folder
  `${CMAKE_INSTALL_PREFIX}\bin` where `CMAKE_INSTALL_PREFIX` is specified
  on the `cmake` command (see above). The binaries will then be installed in
  a subfolder named after the current version stored in the `cfg\VERSION` file.
