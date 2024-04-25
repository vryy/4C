/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for turbulence

\level 2


*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_TURBULENCE_HPP
#define FOUR_C_INPAR_TURBULENCE_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
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

}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
