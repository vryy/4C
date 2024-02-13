/*----------------------------------------------------------------------*/
/*! \file

\brief entry point for cardiac monodomain scalar transport problems

\level 2


*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_CARDIAC_MONODOMAIN_DYN_HPP
#define BACI_SCATRA_CARDIAC_MONODOMAIN_DYN_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

/*! entry point for the solution of electrochemistry problems */
void scatra_cardiac_monodomain_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI cardiac monodomain-module logo on the screen */
void printheartlogo();


BACI_NAMESPACE_CLOSE

#endif  // SCATRA_CARDIAC_MONODOMAIN_DYN_H
