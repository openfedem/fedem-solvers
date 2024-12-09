/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file releaseSingeltons.C
  \brief Cleaning of global heap-allocated objects on program exit.
*/

#include "vpmCommon/SignalHandler.H"
#include "FFlLib/FFlInit.H"
#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"


/*!
  \brief Helper function to release all global heap-allocated objects on exit.
  \param[in] lastPart If \e true, also release the FE library singeltons.

  \details This function is used to deallocate the singelton objects that may
  have been allocated on program exit, mostly to make sure the heap is clean
  when exiting and enable easier detection of other (real) memory leaks.
  Set \a lastPart to \e false if the reducePart() function is going to be
  invoked again in the same session. The FFl and FFaCmdLineArg objects are
  then not released as they are only initialized on startup.
*/

void releaseSingeltons (bool lastPart)
{
  if (lastPart)
  {
    FFl::releaseAllReaders();
    FFl::releaseAllElements();
    FFaCmdLineArg::removeInstance();
    SignalHandler::reset();
  }
  FFaMsg::setMessager();
}
