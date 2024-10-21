// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COUPLING_VOLMORTAR_DEFINES_HPP
#define FOUR_C_COUPLING_VOLMORTAR_DEFINES_HPP

#include "4C_config.hpp"

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
