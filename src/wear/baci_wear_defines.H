/*----------------------------------------------------------------------*/
/*! \file

\brief  norms and tolerances for wear algorithm

\level 1

*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 12/13 |
 *----------------------------------------------------------------------*/
#ifndef BACI_WEAR_DEFINES_H
#define BACI_WEAR_DEFINES_H


/************************************************************************/
/* Wear tolerances                                                      */
/************************************************************************/

// advection map tolerance
#define WEARADVMAP 1.0e-8 /* tolerance for advection map*/

// tolerance for singular matrix
#define WEARSING 1.0e-12 /* singularity check within advection map*/

// convergence check
#define WEARCONV 1.0e-12 /* tolerance for convergence check*/

BACI_NAMESPACE_OPEN

BACI_NAMESPACE_CLOSE

#endif  // WEAR_DEFINES_H
