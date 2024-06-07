/*-----------------------------------------------------------*/
/*! \file

\brief Data container class for geometric search parameters

\level 3

*/
/*-----------------------------------------------------------*/

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

  verbosity_ = Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(io_params, "VERBOSITY");

  write_visualization_ = Core::UTILS::IntegralValue<int>(
      geometric_search_params, "WRITE_GEOMETRIC_SEARCH_VISUALIZATION");
}
FOUR_C_NAMESPACE_CLOSE
