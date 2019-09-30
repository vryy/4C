/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for twophaseflow

\level 2

\maintainer Martin Kronbichler

*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_TWOPHASE_H
#define INPAR_TWOPHASE_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace TWOPHASE
  {
    /// type of surface tension approximation
    enum SurfaceTensionApprox
    {
      surface_tension_approx_none,
      surface_tension_approx_fixed_curvature,
      surface_tension_approx_divgrad_normal,
      surface_tension_approx_nodal_curvature,
      surface_tension_approx_laplacebeltrami
    };

    /// type of surface tension approximation
    enum Transport_Directions
    {
      transport_dir_all,
      transport_dir_normal
    };

    /// type of surface tension approximation
    enum SmoothGradPhi
    {
      // smooth_grad_phi_none,
      smooth_grad_phi_meanvalue,
      smooth_grad_phi_l2_projection,
      smooth_grad_phi_superconvergent_patch_recovery_3D,
      smooth_grad_phi_superconvergent_patch_recovery_2Dz
    };

    enum NodalCurvatureCalc
    {
      l2_projected,
      averaged
    };

    enum LaplaceBeltramiCalc
    {
      matrix_non_smoothed,
      matrix_smoothed,
      matrix_mixed_smoothed
    };

    /// set the twophase parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace TWOPHASE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_TWOPHASE_H
