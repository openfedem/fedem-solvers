// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file testRestart.C

  \brief Test program for restarting the dynamics solver from an in-core state.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 16 Jan 2017
*/

#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>

#include "../solverInterface.h"

int compareResponse (const char*, const char*, double, int = 0);


/*!
  \brief Test program for restarting the dynamics solver from an in-core state.

  \callgraph
*/

int main (int argc, char** argv)
{
  // Lamda function for removing <n> entries from the command-line argument list
  auto&& removeArgs = [&argc,argv](int first, int n)
  {
    argc -= n;
    for (int i = first; i < argc; i++)
      argv[i] = argv[i+n];
  };

  // Check for specification of time window length
  int i, j, nStep = 100;
  std::vector<int> fId;
  const char* vfy = NULL;
  double eps = 0.0;

  for (i = 1; i < argc;)
    if (!strncmp(argv[i],"-Win",4) && i+1 < argc)
    {
      // Extract "-Window <nStep>" from the command-line argument list
      nStep = atoi(argv[i+1]);
      removeArgs(i,2);
    }
    else if (!strcmp(argv[i],"-output"))
    {
      // Extract "-output <f1> <f2> ..." from the command-line argument list
      for (j = i+1; j < argc; j++)
        if (isdigit(argv[j][0]))
          fId.push_back(atoi(argv[j]));
        else
          break;
      removeArgs(i,1+fId.size());
    }
    else if (!strcmp(argv[i],"-verify") && i+2 < argc)
    {
      // Extract "-verify <filename> <eps>" from the command-line argument list
      vfy = argv[i+1];
      eps = atof(argv[i+2]);
      removeArgs(i,3);
    }
    else
      i++;

  // Read input files, preprocess the model and setup the initial configuration
  int status = solverInit(argc,argv);
  if (status) return status; // Simulation failed, aborting...

  // Open file for response outputs
  const int nOut = fId.size();
  std::vector<double> outputs;
  std::ofstream os;
  if (nOut > 0 && vfy)
  {
    outputs.resize(nOut*nStep);
    os.open("outputs.asc");
    os.precision(9);
  }

  int stateSize = getStateSize();
  double* state = new double[stateSize];
  int gagesSize = getGagesSize();
  double* gages = NULL;

  // Caution: We here assume the time step size in the model is constant
  double t  = getCurrentTime();
  double dt = getNextTime(&status) - t;

  // Invoke the dynamics solver for the time windows
  while (solveWindow(nStep,0,0,nOut,fId.data(),
                     NULL,NULL,outputs.data(),0,NULL,&status))
  {
    if (status < 0) return status; // Simulation failed, aborting...

    // Save current state to core array
    if (!saveState(state,stateSize)) return -1;

    if (gagesSize > 0 && !gages)
    {
      // Save the initial gage strains to core array
      gages = new double[gagesSize];
      if (!saveGages(gages,gagesSize)) return -1;
    }

    // Simulation finished, terminate by closing down the result database, etc.
    int dstat = solverDone(false);
    if (status)
      return status;
    else if (dstat)
      return dstat;

    // Print the outputs
    if (!outputs.empty())
      for (i = 1; i <= nStep; i++)
      {
        os << t+dt*i;
        for (j = 0; j < nOut; j++)
          os <<" "<< outputs[nOut*(i-1)+j];
        os << std::endl;
      }

    // Restart the simulation, reading the state and gages from core arrays
    status = solverInit(argc,argv,NULL,state,stateSize,gages,gagesSize);
    if (status) return status; // Simulation failed, aborting...

    if (!outputs.empty())
    {
      t  = getCurrentTime();
      dt = getNextTime(&status) - t;
    }
  }
  delete[] state;
  delete[] gages;

  // Simulation finished, terminate by closing down the result database, etc.
  int dstat = status < 0 ? 0 : solverDone();
  if (status)
    return status;
  else if (dstat)
    return dstat;

  // Verify the simulation by comparing with some reference data
  return compareResponse("outputs.asc",vfy,eps);
}
