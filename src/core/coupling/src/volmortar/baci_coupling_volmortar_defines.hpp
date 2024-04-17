/*----------------------------------------------------------------------*/
/*! \file

\level 1


*----------------------------------------------------------------------*/

#ifndef FOUR_C_COUPLING_VOLMORTAR_DEFINES_HPP
#define FOUR_C_COUPLING_VOLMORTAR_DEFINES_HPP

#include "baci_config.hpp"

FOUR_C_NAMESPACE_OPEN

/************************************************************************/
/* Mortar algorithm parameters                                          */
/************************************************************************/

// MORTAR INTEGRATION
#define VOLMORTARINTTOL 1.0e-12 /* tolerance for assembling gp-values*/

// GEOMETRIC TOLERANCES
#define VOLMORTARELETOL 1.0e-12
#define VOLMORTARCUTTOL 1.0e-12
#define VOLMORTARCUT2TOL -1.0e-12

FOUR_C_NAMESPACE_CLOSE

#endif
