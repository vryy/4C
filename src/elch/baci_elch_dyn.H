/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for Electrochemistry module.

\level 2


*----------------------------------------------------------------------*/

#ifndef BACI_ELCH_DYN_H
#define BACI_ELCH_DYN_H

#include "baci_config.H"

BACI_NAMESPACE_OPEN

/*! entry point for the solution of electrochemistry problems */
void elch_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI electrochemistry-module logo on the screen */
void printlogo();


BACI_NAMESPACE_CLOSE

#endif  // ELCH_DYN_H
