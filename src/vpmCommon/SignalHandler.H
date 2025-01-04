// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file SignalHandler.H

  \brief Handling of various OS signals casing the program to abort.
*/

#ifndef SIGNAL_HANDLER_H
#define SIGNAL_HANDLER_H

#include <string>

//! \brief Function pointer.
typedef void (*VoidFuncPtr)();


/*!
  \brief Static class with methods for signal handling.
*/

class SignalHandler
{
public:
  //! \brief Initializes the signal handler.
  //! \param[in] program Name of executable
  //! \param[in] f1 Function to invoke when exiting without saving
  //! \param[in] f2 Function to invoke when exiting with save
  //! \param[in] f3 Function to invoke when sub-processes should go away also
  static void init (const char* program,
		    VoidFuncPtr f1 = NULL,
		    VoidFuncPtr f2 = NULL,
		    VoidFuncPtr f3 = NULL);
  //! \brief Resets the signal handler to default behaviour.
  static void reset ();

protected:
  //! \brief Invokes the appropriate exit method based on the given signal.
  static void processSignal(int sig);

  //! Function to invoke when exiting without saving
  static VoidFuncPtr emergencyExit;
  //! Function to invoke when exiting with save
  static VoidFuncPtr emergencySaveExit;
  //! Function to invoke when sub-processes should go away also
  static VoidFuncPtr killAllChilds;

  static std::string myProgName; //!< Name of executable
};

#endif
