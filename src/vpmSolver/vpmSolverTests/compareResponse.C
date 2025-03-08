// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file compareResponse.C

  \brief Utility for comparing two ASCII files with numerical data.

  \details This is only used to facilitate regression testing of the solver,
  where a simulation generates some responses (time history data),
  which are to be compared against reference responses stored on some file.
  It is typically invoked at the end of the simulation, after the response to
  be checked has been exported to an ASCII-file.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 13 Dec 2016
*/

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>


/*!
  \brief Utility function to compare two ASCII files with numerical data.
  \param[in] fn1 Path to the first file
  \param[in] fn2 Path to the second file (with reference data)
  \param[in] eps Comparison tolerance
  \param[in] skip Number of initial lines to skip comparison for

  \callergraph
*/

int compareResponse (const char* fn1, const char* fn2, double eps, int skip)
{
  if (!fn1 || !fn2)
    return 0; // Nothing to compare, silence

  std::ifstream f1(fn1);
  if (!f1)
  {
    perror(fn1);
    return -1;
  }

  std::ifstream f2(fn2);
  if (!f2)
  {
    perror(fn2);
    return -2;
  }

  // Lambda function for reading one line of doubles from a stream.
  auto&& readLine = [](std::istream& is, std::vector<double>& values) -> int
  {
    std::string cline;
    if (!std::getline(is,cline))
      return -1; // End-of-file
    else if (cline[0] == '#')
      return 0; // Ignore comment lines

    values.clear();
    size_t itab = 0;
    size_t jtab = cline.find_first_of("\t ,");
    while (jtab < std::string::npos)
    {
      if (itab < jtab)
        values.push_back(atof(cline.substr(itab,jtab-itab).c_str()));
      itab = jtab + 1;
      jtab = cline.find_first_of("\t ,",itab);
    }
    values.push_back(atof(cline.substr(itab).c_str()));
    return values.size();
  };

  // Lambda function for comparing two double values with some tolerance.
  auto&& compareValues = [eps](double v1, double v2) -> bool
  {
    double diff = fabs(v1-v2);
    double refv = fabs(v1+v2)*eps*0.5;
    return diff < refv || diff < eps;
  };

  std::cout <<"\n   * Verifying simulation response in "<< fn1
            <<"\n     against reference data in "<< fn2 << std::endl;
  int numErrors = 0;
  int line, nc1 = 0, nc2 = 0;
  std::vector<double> v1, v2;
  for (line = 1; nc1 >= 0; line++)
    if ((nc1 = readLine(f1,v1)) > 0)
    {
      while ((nc2 = readLine(f2,v2)) == 0);
      if (nc2 < 0)
      {
        std::cerr <<" *** Too few data lines ("<< line
                  <<") in the reference file\n     "<< fn2 << std::endl;
        return -3;
      }
      for (int i = 0; line > skip && i < nc1 && i < nc2; i++)
        if (!compareValues(v1[i],v2[i]))
        {
          numErrors++;
          std::cout <<"  ** Discrepancy in line,column "<< line <<","<< i+1
                    <<": "<< v1[i]-v2[i] <<"\tvalue="<< v1[i]
                    <<"\treference="<< v2[i] << std::endl;
        }
    }

  if (numErrors > 0)
    std::cout <<" *** A total of "<< numErrors
              <<" discrepancies were detected."<< std::endl;
  else
    std::cout <<"   * OK"<< std::endl;

  return numErrors;
}


/*!
  \brief Utility function writing the beginning and end of a file to console.
  \param[in] fileName Path to the file to bewritten
  \param[in] nc1 Number of lines to write from beginning
  \param[in] nc2 Number of lines to write before end

  \callergraph
*/

void writeFile (const char* fileName, int nc1, int nc2)
{
  std::ifstream is(fileName);
  if (!is) return;

  std::cout <<"   * Here is the content of file "<< fileName <<":"<< std::endl;

  int nline = 0;
  std::string cline;
  while (std::getline(is,cline))
    if (++nline <= nc1)
      std::cout << nline <<"\t"<< cline <<"\n";
  if (nline <= nc1) return;

  is.clear();
  is.seekg(0);
  if (nline > nc1+nc2)
    std::cout <<"\n\t...\n\n";
  for (int line = 1; std::getline(is,cline); line++)
    if (line > nc1 && nline-line < nc2)
      std::cout << line <<"\t"<< cline <<"\n";
}
