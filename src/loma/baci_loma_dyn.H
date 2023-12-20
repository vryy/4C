/*---------------------------------------------------------------------*/
/*! \file

\brief Control routine for low-Mach-number flow module.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LOMA_DYN_H
#define BACI_LOMA_DYN_H

#include "baci_config.H"

BACI_NAMESPACE_OPEN

/*! entry point for the solution of low-Mach-number flow problems */
void loma_dyn(int restart /* do we have to perform a restart?  */
);


BACI_NAMESPACE_CLOSE

#endif  // LOMA_DYN_H
