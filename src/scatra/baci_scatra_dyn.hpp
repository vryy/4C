/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport control algorithm
\level 1


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_DYN_HPP
#define FOUR_C_SCATRA_DYN_HPP

#include "baci_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*! entry point for the solution of scalar transport problems */
void scatra_dyn(int restart /* do we have to perform a restart?  */
);

FOUR_C_NAMESPACE_CLOSE

#endif
