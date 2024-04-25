/*----------------------------------------------------------------------*/
/*! \file
 \brief entry point (global control routine) for porous multiphase flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_DYN_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_DYN_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN


/*! entry point for the solution of Lubrication problems */
void porofluidmultiphase_dyn(int restart /* do we have to perform a restart?  */
);


FOUR_C_NAMESPACE_CLOSE

#endif
