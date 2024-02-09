/*---------------------------------------------------------------------*/
/*! \file

\brief Control routine for low-Mach-number flow module.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LOMA_DYN_HPP
#define BACI_LOMA_DYN_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

/*! entry point for the solution of low-Mach-number flow problems */
void loma_dyn(int restart /* do we have to perform a restart?  */
);


BACI_NAMESPACE_CLOSE

#endif  // LOMA_DYN_H
