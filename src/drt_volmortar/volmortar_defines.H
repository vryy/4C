/*!----------------------------------------------------------------------

\level 1

\maintainer Alexander Popp

*----------------------------------------------------------------------*/

#ifndef VOLMORTAR_DEFINES_H_
#define VOLMORTAR_DEFINES_H_

/************************************************************************/
/* Mortar algorithm parameters                                          */
/************************************************************************/

// MORTAR INTEGRATION
#define VOLMORTARINTTOL 1.0e-12 /* tolerance for assembling gp-values*/

// GEOMETRIC TOLERANCES
#define VOLMORTARELETOL 1.0e-12
#define VOLMORTARCUTTOL 1.0e-12
#define VOLMORTARCUT2TOL -1.0e-12

#endif /* VOLMORTAR_DEFINES_H_ */
