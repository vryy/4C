/*----------------------------------------------------------------------*/
/*! \file

\level 1


*----------------------------------------------------------------------*/

#ifndef COUPLING_VOLMORTAR_DEFINES_H
#define COUPLING_VOLMORTAR_DEFINES_H

/************************************************************************/
/* Mortar algorithm parameters                                          */
/************************************************************************/

// MORTAR INTEGRATION
#define VOLMORTARINTTOL 1.0e-12 /* tolerance for assembling gp-values*/

// GEOMETRIC TOLERANCES
#define VOLMORTARELETOL 1.0e-12
#define VOLMORTARCUTTOL 1.0e-12
#define VOLMORTARCUT2TOL -1.0e-12

#endif  // COUPLING_VOLMORTAR_DEFINES_H
