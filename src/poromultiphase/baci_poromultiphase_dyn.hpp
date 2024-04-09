/*----------------------------------------------------------------------*/
/*! \file
 \brief entry point (global control routine) for poroelasticity multiphase flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_DYN_HPP
#define FOUR_C_POROMULTIPHASE_DYN_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

/*! entry point for the solution of poro multiphase problems problems */
void poromultiphase_dyn(int restart /* do we have to perform a restart?  */
);

BACI_NAMESPACE_CLOSE

#endif
