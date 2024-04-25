/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for Electrochemistry module.

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_ELCH_DYN_HPP
#define FOUR_C_ELCH_DYN_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*! entry point for the solution of electrochemistry problems */
void elch_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI electrochemistry-module logo on the screen */
void printlogo();


FOUR_C_NAMESPACE_CLOSE

#endif
