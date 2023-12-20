/*--------------------------------------------------------------------------*/
/*! \file

\brief entry point for the solution of Lubrication problems

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef BACI_LUBRICATION_DYN_H
#define BACI_LUBRICATION_DYN_H

#include "baci_config.H"

BACI_NAMESPACE_OPEN


/*! entry point for the solution of Lubrication problems */
void lubrication_dyn(int restart /* do we have to perform a restart?  */
);


BACI_NAMESPACE_CLOSE

#endif  // LUBRICATION_DYN_H
