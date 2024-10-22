// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_DEFINES_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_DEFINES_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/************************************************************************/
/* Projection                                                           */
/************************************************************************/
#define PROJMAXITER 10 /* max number of steps for local Newton */
// default: 10
#define PROJOUTPUT false /* output of projection algorithm */
// default: false
#define XIETATOL                                                                               \
  1.0e-9 /* used for check if two coordinates xi, eta in element parameter space are identical \
default: 1.0e-9  --> segments (in parameter space) cannot be smaller than this */
#define COLINEARTOL 1.0e-08 /* used for check if elements are colinear */
// default: 1.0e-08
#define CONVTOLNEWTONPROJ 1.0e-10 /* tolerance for Newton projection algorithm */
// default: 1.0e-10
/************************************************************************/
/* Random                                                               */
/************************************************************************/
#define KAPPAINVTOL 1.0e-10 /* used if kappa vector can be inverted */
// default: 1.0e-10
#define CONVTOLNEWTONMESHMOVEMENT                                                                                       \
  1.0e-09                      /* used for convergence check of Newton for                                              \
//default: 1.0e-9                             recomputation of eta and xi in current configuration \
                     */
#define MESHMOVEMENTMAXITER 10 /* max number of steps for local Newton (mesh movement) */

FOUR_C_NAMESPACE_CLOSE

#endif
