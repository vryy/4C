/*-----------------------------------------------------------*/
/*! \file

\brief Data container class for geometric search parameters

\level 3

*/
/*-----------------------------------------------------------*/

#include <Teuchos_ParameterList.hpp>

#include "baci_lib_globalproblem.H"

#include "baci_discretization_geometric_search_params.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::GEOMETRICSEARCH::GeometricSearchParams::GeometricSearchParams()
    : beam_radius_extension_factor_(-1), sphere_radius_extension_factor_(-1)
{
  Teuchos::ParameterList const& params_list = ::DRT::Problem::Instance()->GeometricSearchParams();

  beam_radius_extension_factor_ = params_list.get<double>("BEAM_RADIUS_EXTENSION_FACTOR");
  dsassert(!std::signbit(beam_radius_extension_factor_),
      "Beam radius extension factor needs to be positive!");

  sphere_radius_extension_factor_ = params_list.get<double>("SPHERE_RADIUS_EXTENSION_FACTOR");
  dsassert(!std::signbit(sphere_radius_extension_factor_),
      "Sphere radius extension factor needs to be positive!");

  Teuchos::ParameterList const& params_list_io = ::DRT::Problem::Instance()->IOParams();
  verbosity_ = ::DRT::INPUT::IntegralValue<IO::verbositylevel>(params_list_io, "VERBOSITY");
}