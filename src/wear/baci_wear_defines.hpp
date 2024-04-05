/*----------------------------------------------------------------------*/
/*! \file

\brief  norms and tolerances for wear algorithm

\level 1

*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 12/13 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_WEAR_DEFINES_HPP
#define FOUR_C_WEAR_DEFINES_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

/************************************************************************/
/* Wear tolerances                                                      */
/************************************************************************/

// advection map tolerance
#define WEARADVMAP 1.0e-8 /* tolerance for advection map*/

// tolerance for singular matrix
#define WEARSING 1.0e-12 /* singularity check within advection map*/

// convergence check
#define WEARCONV 1.0e-12 /* tolerance for convergence check*/

BACI_NAMESPACE_CLOSE

#endif  // WEAR_DEFINES_H
