/*---------------------------------------------------------------------*/
/*! \file

\brief Control routine for low-Mach-number flow module.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_LOMA_DYN_HPP
#define FOUR_C_LOMA_DYN_HPP

#include "baci_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*! entry point for the solution of low-Mach-number flow problems */
void loma_dyn(int restart /* do we have to perform a restart?  */
);


FOUR_C_NAMESPACE_CLOSE

#endif
