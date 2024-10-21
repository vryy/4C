// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_geometric_search_params.hpp"

#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::GeometricSearch::GeometricSearchParams::GeometricSearchParams(
    const Teuchos::ParameterList& geometric_search_params, const Teuchos::ParameterList& io_params)
    : beam_radius_extension_factor_(-1), sphere_radius_extension_factor_(-1)
{
  beam_radius_extension_factor_ =
      geometric_search_params.get<double>("BEAM_RADIUS_EXTENSION_FACTOR");
  FOUR_C_ASSERT(!std::signbit(beam_radius_extension_factor_),
      "Beam radius extension factor needs to be positive!");

  sphere_radius_extension_factor_ =
      geometric_search_params.get<double>("SPHERE_RADIUS_EXTENSION_FACTOR");
  FOUR_C_ASSERT(!std::signbit(sphere_radius_extension_factor_),
      "Sphere radius extension factor needs to be positive!");

  verbosity_ = Teuchos::getIntegralValue<Core::IO::Verbositylevel>(io_params, "VERBOSITY");

  write_visualization_ = geometric_search_params.get<bool>("WRITE_GEOMETRIC_SEARCH_VISUALIZATION");
}
FOUR_C_NAMESPACE_CLOSE
