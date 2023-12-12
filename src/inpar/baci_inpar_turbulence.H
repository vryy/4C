/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for turbulence

\level 2


*/

/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_TURBULENCE_H
#define BACI_INPAR_TURBULENCE_H

#include "baci_inpar_parameterlist_utils.H"

BACI_NAMESPACE_OPEN

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
BACI_NAMESPACE_CLOSE

#endif  // INPAR_TURBULENCE_H
