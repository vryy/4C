// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_TURBULENCE_HPP
#define FOUR_C_INPAR_TURBULENCE_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace FLUID
  {
    //------------------------------------------------------------------
    // parameters for multifractal subgrid-scale modeling
    //------------------------------------------------------------------
    //! scale separation
    enum ScaleSeparation
    {
      no_scale_sep,
      box_filter,
      algebraic_multigrid_operator
    };

    //! definition Reynolds number: reference length
    enum RefLength
    {
      cube_edge,
      sphere_diameter,
      streamlength,
      gradient_based,
      metric_tensor
    };

    //! definition Reynolds number: reference velocity
    enum RefVelocity
    {
      strainrate,
      resolved,
      fine_scale
    };

  }  // namespace FLUID

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
