/*----------------------------------------------------------------------*/
/*! \file

\brief entry point for cardiac monodomain scalar transport problems

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_CARDIAC_MONODOMAIN_DYN_HPP
#define FOUR_C_SCATRA_CARDIAC_MONODOMAIN_DYN_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*! entry point for the solution of electrochemistry problems */
void scatra_cardiac_monodomain_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI cardiac monodomain-module logo on the screen */
void printheartlogo();


FOUR_C_NAMESPACE_CLOSE

#endif
