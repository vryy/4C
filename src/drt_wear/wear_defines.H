/*----------------------------------------------------------------------*/
/*!
\file wear_defines.H

\brief  norms and tolerances for wear algorithm

\level 1

<pre>

\maintainer Alexander Popp

</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 12/13 |
 *----------------------------------------------------------------------*/
#ifndef WEAR_DEFINES_H_
#define WEAR_DEFINES_H_


/************************************************************************/
/* Wear tolerances                                                      */
/************************************************************************/

// advection map tolerance
#define WEARADVMAP 1.0e-8 /* tolerance for advection map*/

// tolerance for singular matrix
#define WEARSING 1.0e-12 /* singularity check within advection map*/

// convergence check
#define WEARCONV 1.0e-12 /* tolerance for convergence check*/

#endif /* WEAR_DEFINES_H_ */
