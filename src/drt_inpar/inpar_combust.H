/*----------------------------------------------------------------------*/
/*!
\file inpar_combust.H

\brief Input parameters for combustion

\level 3

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_COMBUST_H
#define INPAR_COMBUST_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  namespace COMBUST
  {
    /// type of combustion problem
    enum CombustionType
    {
      combusttype_premixedcombustion,
      combusttype_twophaseflow,
      combusttype_twophaseflow_surf,
      combusttype_twophaseflowjump
    };

    /// type of reinitialization for level set function (G-function)
    enum ReInitialActionGfunc
    {
      reinitaction_none,
      reinitaction_byfunction,
      reinitaction_signeddistancefunction,
      reinitaction_fastsigneddistancefunction,
      reinitaction_sussman
    };

    /// type of domain and boundary integration cell construction
    enum XFEMIntegration
    {
      xfemintegration_cut,
      xfemintegration_tetrahedra,
      xfemintegration_hexahedra
    };

    /// type of xfem time integration
    enum XFEMTimeIntegration
    {
      xfemtimeint_donothing,
      xfemtimeint_semilagrange,
      xfemtimeint_extrapolationold,
      xfemtimeint_extrapolationnew,
      xfemtimeint_mixedSLExtrapol,
      xfemtimeint_mixedSLExtrapolNew,
      xfemtimeint_mixedghostSL,
      xfemtimeint_mixedghostExtrapol,
      xfemtimeint_mixedghostSLExtrapol
    };

    /// type of xfem time integration of enrichment values
    enum XFEMTimeIntegrationEnr
    {
      xfemtimeintenr_donothing,
      xfemtimeintenr_quasistatic,
      xfemtimeintenr_project,
      xfemtimeintenr_project_scalar,
      xfemtimeintenr_extrapolate,
      xfemtimeintenr_extrapolate_scalar
    };

    /// type of xfem time integration of enrichment values computation
    enum XFEMTimeIntegrationEnrComp
    {
      xfemtimeintenr_full,
      xfemtimeintenr_standard,
      xfemtimeintenr_minimal
    };

    /// enriched fields when using kink enrichments
    enum SelectedEnrichment
        {
          selectedenrichment_both,
          selectedenrichment_velocity,
          selectedenrichment_pressure,
          selectedenrichment_none
        };

    /// type of initial fluid field
    enum InitialField
    {
      initfield_zero_field,
      initfield_field_by_function,
      initfield_disturbed_field_by_function,
      initfield_flame_vortex_interaction,
      initfield_darrieus_landau_instability,
      initfield_beltrami_flow
    };

    /// type of surface tension approximation
    enum SurfaceTensionApprox
    {
      surface_tension_approx_none,
      surface_tension_approx_fixed_curvature,
      surface_tension_approx_divgrad,
      surface_tension_approx_divgrad_normal,
      surface_tension_approx_nodal_curvature,
      surface_tension_approx_laplacebeltrami,
      surface_tension_approx_laplacebeltrami_smoothed
    };

    /// type of surface tension approximation
    enum SmoothGradPhi
    {
      smooth_grad_phi_none,
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

    /// type of velocity jump
    enum VelocityJumpType
    {
      vel_jump_none,
      vel_jump_const,
      vel_jump_premixed_combustion
    };

    /// type of velocity jump
    enum FluxJumpType
    {
      flux_jump_none,
      flux_jump_const,
      flux_jump_premixed_combustion,
      flux_jump_surface_tension
    };

    /// analyical solution for Nitsche
    enum NitscheError
    {
      nitsche_error_none,
      nitsche_error_static_bubble_nxnx1,
      nitsche_error_static_bubble_nxnxn,
      nitsche_error_shear,
      nitsche_error_couette_20x20x1,
      nitsche_error_straight_bodyforce,
      nitsche_error_ellipsoid_bubble_2D,
      nitsche_error_ellipsoid_bubble_3D,
      nitsche_error_beltrami
    };

    /// type of velocity jump
    enum WeightType
    {
      weight_visc_based_harmonic,
      weight_intersection_visc_based_harmonic
    };

    /// set the combust parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace COMBUST

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_COMBUST_H

