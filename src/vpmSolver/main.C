/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file vpmSolver/main.C

  \brief This file contains the C++ main program for the FEDEM dynamics solver.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 2 Dec 2016
*/

#include "solverInterface.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

int compareResponse (const char*, const char*, double, int);
void writeFile (const char*, int, int);


/*!
  \brief Main program for the FEDEM dynamics solver.

  \details The main program contains very little logic.
  It uses functions from the solverInterface API to initialize and solve the
  dynamic problem at each time step. It can also invoke a response verification
  before program termination, in which the calculated response is compared with
  some reference data. This is mainly used to set up regression tests.

  \callgraph
*/

int main (int argc, char** argv)
{
  // Check if a file for response verification was specified
  const char* response = NULL;
  const char* verify = NULL;
  double epsTol = 0.0;
  int skipIL = 0;
  if (argc > 3 && !strcmp(argv[argc-3],"-verify"))
  {
    verify = argv[argc-2];
    epsTol = atof(argv[argc-1]);
    argc -= 3;
  }
  else if (argc > 4 && !strcmp(argv[argc-4],"-verify"))
  {
    verify = argv[argc-3];
    epsTol = atof(argv[argc-2]);
    skipIL = atoi(argv[argc-1]);
    argc -= 4;
  }

  char resfile[128] = "\0";
  char crvfile[128] = "\0";

  // Lambda function printing a console error message on failure.
  // It also prints portions of the fedem_solver.res file.
  auto&& failure = [argv,verify,resfile](const char* prg, int stat)
  {
    std::cerr <<" *** "<< argv[0] <<": "<< prg
              <<" failed ("<< stat <<")"<< std::endl;
    if (verify) writeFile(resfile,200,100);
    return stat;
  };

  // Read input files, preprocess the model and set up the initial configuration
  int status = solverInit(argc,argv);
  if (status) return failure("solverInit",status);

  // Get path to the fedem_solver.res file for this run
  getFileName("resfile",resfile,128);
  if (verify) // Get path to the exported curves file (if any) for this run
    response = getFileName("curvePlotFile",crvfile,128);

  // Time step loop.
  // Invoke the solver step-by-step until specified end time is reached,
  // or an error occur.
  while (solveNext(&status))
    if (status < 0)
      return failure("solveNext",status); // Simulation failed, aborting...

  // Simulation finished, terminate by closing down the result database, etc.
  int dstat = solverDone();
  if (status)
    return failure("solveNext",status);
  else if (dstat)
    return failure("solverDone",dstat);

  // Verify exported curve data against the reference data, if specified
  if ((status = compareResponse(response,verify,epsTol,skipIL)))
  {
    std::cerr <<" *** "<< argv[0]
              <<": Comparison failed ("<< status <<")"<< std::endl;
    writeFile(resfile,200,100);
  }

  return status;
}
