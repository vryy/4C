/*----------------------------------------------------------------------*/
/*! \file
\brief entry point for level-set transport problems
\level 2

*----------------------------------------------------------------------*/

#ifndef FOUR_C_LEVELSET_DYN_HPP
#define FOUR_C_LEVELSET_DYN_HPP

#include "baci_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*! entry point for the solution of level-set problems */
void levelset_dyn(int restart /* do we have to perform a restart?  */
);

FOUR_C_NAMESPACE_CLOSE

#endif
