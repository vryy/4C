/*----------------------------------------------------------------------*/
/*!
\brief Input parameters for turbulence

\level 2

\maintainer Martin Kronbichler

*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_TURBULENCE_H
#define INPAR_TURBULENCE_H

#include "inpar_parameterlist_utils.H"

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
#endif  // #ifndef INPAR_TURBULENCE_H
