/*---------------------------------------------------------------------*/
/*! \file
\brief Control routine for structure with ale problems (wear).

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_WEAR_DYN_H
#define BACI_WEAR_DYN_H

#include "baci_config.H"

BACI_NAMESPACE_OPEN

/* entry point for the solution of structural ale wear problems */
void wear_dyn_drt(int restart);

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // WEAR_DYN_H
