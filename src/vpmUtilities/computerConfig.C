/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file computerConfig.C
  \brief Global functions for extracting information on the running computer.
  \details This file contains the implementation of the following functions
  that can be invoked from Fortran programs:

  - computerconfiginterface::getComputerConfig
  - computerconfiginterface::getUserName

  No further documentation is provided here.
  The methods are documented in the computerConfigInterface.f90 file.
*/

#include "FFaLib/FFaOS/FFaFortran.H"

#if defined(win32) || defined(win64)
#include <windows.h>
#else
#include <sys/utsname.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>


namespace
{
  // \brief Helper to append a string to a fixed-size character buffer.
  void appends (char* dest, const char* src, const int nchar)
  {
    if (int mchar = nchar - strlen(dest) - 2; mchar > 0)
      strncat(strcat(dest," "),src,mchar);
  }
}


SUBROUTINE(getcomputerconfig,GETCOMPUTERCONFIG) (char* cid, const int nchar)
{
  size_t l;
#if defined(win32) || defined(win64)
  DWORD cchBuff = BUFSIZ;
  TCHAR tchBuffer[BUFSIZ];
  LPTSTR hn = tchBuffer;
  char* os; char* id;
  GetComputerName(hn,&cchBuff);
  strncpy(cid,hn,nchar);
  os = getenv("OS");
  appends(cid,os,nchar);
  id = getenv("PROCESSOR_IDENTIFIER");
  appends(cid,id,nchar);
#else
  struct utsname name;
  uname(&name);
  strncpy(cid,name.nodename,nchar);
  appends(cid,name.sysname,nchar);
  appends(cid,name.release,nchar);
  appends(cid,name.version,nchar);
  appends(cid,name.machine,nchar);
#endif
  l = strlen(cid);
  if ((int)l < nchar) memset(cid+l,' ',nchar-l);
}


SUBROUTINE(getusername,GETUSERNAME) (char* cuser, const int nchar)
{
  size_t l;
#if defined(win32) || defined(win64)
  DWORD cchBuff = BUFSIZ;
  TCHAR tchBuffer[BUFSIZ];
  LPTSTR usr = tchBuffer;
  GetUserName(usr,&cchBuff);
#else
  char* usr = getenv("USER");
#endif
  strncpy(cuser, usr ? usr : "(none)", nchar-1);
  l = usr ? strlen(usr) : 6;
  if ((int)l < nchar) memset(cuser+l,' ',nchar-l);
}
