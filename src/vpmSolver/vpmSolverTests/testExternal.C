// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file testExternal.C

  \brief Test program for external functions in the dynamics solver.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 5 Dec 2016
*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "../solverInterface.h"

bool readFsiFile (const char*, const char*, std::string&);
int compareResponse (const char*, const char*, double, int = 0);


/*!
  \brief Hard-coded external function evaluations to test against.
*/

static double extFunc (int id, double x)
{
  switch (id) {
  case 1: return 1.0e7*sin(x);
  case 2: return 5.0e3*(cos(2.5*x)-1.0);
  }
  return 0.0;
}


/*!
  \brief Helper to print out response data to the given stream for a time step.
*/

static void printOut (std::ostream& os, bool formatted,
                      double time, const double* outputs, size_t nOut)
{
  if (nOut < 1) return;
  if (formatted) os <<"Time: ";
  os << time;
  for (size_t j = 0; j < nOut; j++)
    if (formatted)
      os <<", F"<< 1+j <<" = "<< outputs[j];
    else
      os <<" "<< outputs[j];
  os << std::endl;
}


/*!
  \brief Test program for external functions in the dynamics solver.

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

  // Check for testing to pass the fsi-file content in a character variable,
  // if input and/or output functions are to be used, and
  // if a reference file for comparing the response is specified
  int i, j, nIn = 0;
  std::vector<int> fId;
  const char* cwd = NULL;
  const char* fsi = NULL;
  const char* vfy = NULL;
  const char* exf = NULL;
  double eps = 0.0;
  bool tWind = false;

  for (i = 1; i < argc;)
    if (!strcmp(argv[i],"-cwd") && i+1 < argc)
      cwd = argv[++i];

    else if (!strcmp(argv[i],"-fsi") && i+1 < argc)
    {
      // Extract "-fsi <filename>" from the command-line argument list
      fsi = argv[i+1];
      removeArgs(i,2);
    }
    else if (!strcmp(argv[i],"-input") && i+1 < argc)
    {
      // Extract "-input <nIn>" from the command-line argument list
      nIn = atoi(argv[i+1]);
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
    else if (!strcmp(argv[i],"-Window"))
    {
      // Testing the time windows functionality
      tWind = true;
      removeArgs(i,1);
    }
    else if (!strcmp(argv[i],"-writeExtFunc") && i+1 < argc)
    {
      // Write out the external function values to separate file
      exf = argv[i+1];
      removeArgs(i,2);
    }
    else // other arguments to be processed by solverInit()
      i++;

  // Read the whole solver input file into a text string
  std::string chfsi;
  if (!readFsiFile(cwd,fsi,chfsi))
    return -99;

  // Read input files, preprocess the model and setup the initial configuration
  int status = solverInit(argc, argv, chfsi.empty() ? NULL : chfsi.c_str());
  if (status < 0) return status; // Simulation failed, aborting...

  // Extract and print the outputs
  std::vector<double> outputs;
  std::ofstream os, extOs;
  if (fId.size() > 3 || vfy)
  {
    os.open("outputs.asc");
    outputs.reserve(fId.size());
    for (int funcId : fId)
      outputs.push_back(getFunc(funcId));
    if (vfy) os.precision(9);
    printOut(os,false,getCurrentTime(),outputs.data(),fId.size());
  }

  // Write external function values to file
  bool writeCSV = false;
  if (exf && !tWind && nIn >= 2)
  {
    writeCSV = strstr(exf,".csv");
    extOs.open(exf);
    extOs.precision(16);
    if (writeCSV)
      extOs <<"_t,F1,F2\n";
    else
      extOs <<"# External function values test file\n"
            <<"#DESCRIPTION\tF1\tDummy\tF2\tJalla\n";
  }
  else
    exf = NULL;

  // Time step loop
  bool doContinue = true;
  status = 0;
  if (tWind)
  {
    const int nOut = fId.size();
    const int nStep = 100;
    outputs.resize(nOut*nStep);
    // Caution: We here assume the time step size in the model is constant
    double t  = getCurrentTime();
    double dt = getNextTime(&status) - t;

    // Generate inputs
    std::vector<double> inputs;
    for (i = 1; i <= nStep; i++)
      for (j = 1; j <= nIn; j++)
        inputs.push_back(extFunc(j,t+dt*i));

    // Invoke the dynamics solver for the time windows
    solveWindow(nStep,0,nIn,nOut,fId.data(),
                NULL,inputs.data(),outputs.data(),0,NULL,&status);
    if (status < 0) return status; // Simulation failed, aborting...

    // Print the outputs
    for (i = 1; i <= nStep; i++)
      if (os)
	printOut(os,false,t+dt*i,outputs.data()+nOut*(i-1),nOut);
      else
	printOut(std::cout,true,t+dt*i,outputs.data()+nOut*(i-1),nOut);
  }
  else while (doContinue)
  {
    // Evaluate the external functions
    double t = getNextTime(&status);
    for (j = 1; j <= nIn; j++)
      if (!setExtFunc(j,extFunc(j,t)))
        --status;

    long int ms = t*1.0e6+0.1;
    if (writeCSV) // Write external function values to file
      extOs << ms <<","<< extFunc(1,t) <<","<< extFunc(2,t) <<"\n";
    else if (exf) // Write external function values to file
      extOs << t <<"\t"<< extFunc(1,t) <<"\t2.0\t"<< extFunc(2,t) <<"\t4.0\n";

    // Invoke the solver to advance the time one step forward
    if (status == 0)
      doContinue = solveNext(&status) && status >= 0;
    if (status < 0) return status; // Simulation failed, aborting...

    // Extract and print the outputs
    outputs.clear();
    for (int funcId : fId)
      outputs.push_back(getFunc(funcId));
    if (os)
      printOut(os,false,t,outputs.data(),fId.size());
    else
      printOut(std::cout,true,t,outputs.data(),fId.size());
  }

  if (os) os.close();
  if (extOs) extOs.close();

  // Simulation finished, terminate by closing down the result database, etc.
  int dstat = solverDone();
  if (status)
    return status;
  else if (dstat)
    return dstat;

  // Verify the simulation by comparing with some reference data
  return compareResponse("outputs.asc",vfy,eps);
}
