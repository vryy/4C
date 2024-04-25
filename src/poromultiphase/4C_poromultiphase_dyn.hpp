/*----------------------------------------------------------------------*/
/*! \file
 \brief entry point (global control routine) for poroelasticity multiphase flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_DYN_HPP
#define FOUR_C_POROMULTIPHASE_DYN_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*! entry point for the solution of poro multiphase problems problems */
void poromultiphase_dyn(int restart /* do we have to perform a restart?  */
);

FOUR_C_NAMESPACE_CLOSE

#endif
