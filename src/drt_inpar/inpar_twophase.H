/*----------------------------------------------------------------------*/
/*!
\file inpar_twophase.H

\brief Input parameters for combustion

<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
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
      //surface_tension_approx_fixed_curvature,
      //surface_tension_approx_divgrad,
      surface_tension_approx_divgrad_normal,
      surface_tension_approx_nodal_curvature,
      surface_tension_approx_laplacebeltrami
    };

    /// type of surface tension approximation
    enum SmoothGradPhi
    {
      //smooth_grad_phi_none,
      smooth_grad_phi_meanvalue,
      smooth_grad_phi_leastsquares_3D,
      smooth_grad_phi_leastsquares_2Dx,
      smooth_grad_phi_leastsquares_2Dy,
      smooth_grad_phi_leastsquares_2Dz,
      smooth_grad_phi_l2_projection
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

//    /// type of velocity jump
//    enum WeightType
//    {
//      weight_visc_based_harmonic,
//      weight_intersection_visc_based_harmonic
//    };

  }  // namespace TWOPHASE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_TWOPHASE_H

