// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file testInverse.C

  \brief Test program for components of the inverse method.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 23 Jan 2018
*/

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include "../solverInterface.h"

#ifndef M_PI
//! The value of &pi;
#define M_PI 3.14159265358979323846
#endif

int compareResponse (const char*, const char*, double, int = 0);

using Ipair  = std::pair<int,int>; //!< Convenience type definition
using Ipairs = std::vector<Ipair>; //!< Convenience type definition


/*!
  \brief Test program for components of the inverse method.

  \callgraph
*/

int main (int argc, char** argv)
{
  // Lamda function for parsing an int from the command-line argument list
  auto&& parseInt = [argv](int iarg, std::vector<int>& ival)
  {
    if (!isdigit(argv[iarg][0]))
      return false;

    ival.push_back(atoi(argv[iarg]));
    return true;
  };

  // Lamda function for parsing an int pair from the command-line argument list
  auto&& parseInts = [argv](int iarg, Ipairs& ip)
  {
    if (!isdigit(argv[iarg][0]) || !isdigit(argv[iarg+1][0]))
      return false;

    ip.push_back({ atoi(argv[iarg]), atoi(argv[iarg+1]) });
    return true;
  };

  // Lamda function for removing <n> entries from the command-line argument list
  auto&& removeArgs = [&argc,argv](int first, int n)
  {
    argc -= n;
    for (int i = first; i < argc; i++)
      argv[i] = argv[i+n];
  };

  // Lamda function for converting a baseId,DOF pair to an equation number
  auto&& addEquation = [](const Ipair& dof, std::vector<int>& meqn)
  {
    int eqs[6];
    getEquations(dof.first,eqs);
    meqn.push_back(eqs[dof.second-1]);
    std::cout <<"baseID,ldof="<< dof.first <<","<< dof.second
              <<" --> eqn="<< meqn.back() << std::endl;
  };

  // Check for specification of object base ID and local DOF pairs in which the
  // displacements should be measured, and applied loads should be calculated
  // by means of inverse method
  std::vector<int> fId, dId;
  const char* vfy = NULL;
  double eps = 0.0;
  Ipairs d, f;
  int i, j;

  for (i = 1; i < argc;)
    if (!strncmp(argv[i],"-dis",4))
    {
      // Extract "-dis <t1> <d1> <t2> <d2> ..." from the command-line arguments
      for (j = i+1; j+1 < argc && parseInts(j,d); j += 2);
      removeArgs(i,1+2*d.size());
    }
    else if (!strncmp(argv[i],"-forc",5))
    {
      // Extract "-forc <t1> <d1> <t2> <d2> ..." from the command-line arguments
      for (j = i+1; j+1 < argc && parseInts(j,f); j += 2);
      removeArgs(i,1+2*f.size());
    }
    else if (!strcmp(argv[i],"-dfunc"))
    {
      // Extract "-dfunc <f1> <f2> ..." from the command-line argument list
      for (j = i+1; j < argc && parseInt(j,dId); j++);
      removeArgs(i,1+dId.size());
    }
    else if (!strcmp(argv[i],"-output"))
    {
      // Extract "-output <f1> <f2> ..." from the command-line argument list
      for (j = i+1; j < argc && parseInt(j,fId); j++);
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
  if (status != 0) return status; // Simulation failed, aborting...

  // Get the equation numbers corresponding to the specified input
  std::vector<int> meqn_d, meqn_f;
  for (const Ipair& dof : d) addEquation(dof,meqn_d);
  for (const Ipair& dof : f) addEquation(dof,meqn_f);

  // Extract and print the outputs
  std::ofstream os;
  if (!fId.empty() && vfy)
  {
    os.open("outputs.asc");
    os.precision(9);
    os << getCurrentTime();
    for (int funcId : fId) os <<" "<< evalFunc(funcId);
    os << std::endl;
  }
  else
    vfy = NULL;

  // Time step loop
  std::vector<double> dis(d.size(),0.0);
  for (bool doContinue = true; doContinue && status == 0;)
  {
    // Evaluate the displacements (emulated sensor data)
    double t = getNextTime(&status);
    dis[0] = 0.020*sin(2.0*M_PI*t);
    dis[1] = 0.003*sin(2.0*M_PI*t);
    for (size_t k = 0; k < dId.size() && k < dis.size(); k++)
      dis[k] = evalFunc(dId[k],NULL,t);
    std::cout <<"   * Solving at t="<< t <<" for displacements";
    for (double v : dis) std::cout <<" "<< v;
    std::cout << std::endl;

    // Perform the system iterations
    doContinue = solveInverse(dis.data(),meqn_d.data(),meqn_f.data(),
                              meqn_d.size(),meqn_f.size(),&status);
    if (status < 0) return status; // Simulation failed, aborting...

    if (vfy) // Extract and print the outputs
    {
      os << t;
      for (int funcId : fId) os <<" "<< evalFunc(funcId);
      os << std::endl;
    }
  }

  // Simulation finished, terminate by closing down the result database, etc.
  int dstat = solverDone();
  if (status)
    return status;
  else if (dstat)
    return dstat;

  // Verify the simulation by comparing with some reference data
  return compareResponse("outputs.asc",vfy,eps);
}
