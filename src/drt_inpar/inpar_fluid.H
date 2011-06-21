/*----------------------------------------------------------------------*/
/*!
\file inpar_combust.H

\brief Input parameters for fluid

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_FLUID_H
#define INPAR_FLUID_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  namespace FLUID
  {
    //! physical type of the fluid flow (incompressible, Boussinesq approximation, varying density, loma)
    enum PhysicalType
    {
      incompressible,
      boussinesq,
      varying_density,
      loma
    };
    //! flag to (de)activate PSPG stabilization
    enum PSPG
    {
      pstab_assume_inf_sup_stable,  // no pspg stabilization
      pstab_use_pspg                // use pspg stabilization
    };
    //! flag to (de)activate SUPG stabilization
    enum SUPG
    {
      convective_stab_none,
      convective_stab_supg
    };
    //! flag to select the type of viscosity stabilization
    enum VStab
    {
      viscous_stab_none,
      viscous_stab_gls,
      viscous_stab_gls_only_rhs,
      viscous_stab_usfem,
      viscous_stab_usfem_only_rhs
    };
    //! flag to select the type of reactive stabilization
    enum RStab
    {
      reactive_stab_none,
      reactive_stab_gls,
      reactive_stab_usfem
    };
    //! flag to (de)activate continuity stabilization
    enum CStab
    {
      continuity_stab_none,
      continuity_stab_yes
    };
    //! flag to select the type of cross stress stabilization
    enum CrossStress
    {
      cross_stress_stab_none,
      cross_stress_stab,
      cross_stress_stab_only_rhs
    };
    //! flag to select the type of Reynolds stress stabilization
    enum ReynoldsStress
    {
      reynolds_stress_stab_none,
      reynolds_stress_stab,
      reynolds_stress_stab_only_rhs
    };

    enum SubscalesTD
    {
      subscales_none,
      subscales_quasistatic,
      subscales_time_dependent
    };

    enum Transient
    {
      inertia_stab_drop,
      inertia_stab_keep,
      inertia_stab_keep_complete
    };

    enum TauType
    {
      tau_taylor_hughes_zarins,
      tau_taylor_hughes_zarins_wo_dt,
      tau_taylor_hughes_zarins_whiting_jansen,
      tau_taylor_hughes_zarins_whiting_jansen_wo_dt,
      tau_taylor_hughes_zarins_scaled,
      tau_taylor_hughes_zarins_scaled_wo_dt,
      tau_franca_barrenechea_valentin_frey_wall,
      tau_franca_barrenechea_valentin_frey_wall_wo_dt,
      tau_shakib_hughes_codina,
      tau_shakib_hughes_codina_wo_dt,
      tau_codina,
      tau_codina_wo_dt,
      tau_franca_madureira_valentin_badia_codina,
      tau_franca_madureira_valentin_badia_codina_wo_dt,
      tau_not_defined
    };

    enum TauType_genalpha
    {
    smoothed_franca_barrenechea_valentin_wall,
    franca_barrenechea_valentin_wall         ,
    franca_barrenechea_valentin_codina       ,
    fbvw_gradient_based_hk                   ,
    codina                                   ,
    taylor_hughes_zarins_whiting_jansen      ,
    fbvw_wo_dt                               ,
    tautype_not_defined
  };

    /// time integration scheme
    enum TimeIntegrationScheme
    {
      timeint_stationary,
      timeint_one_step_theta,
      timeint_npgenalpha, // fluid3's implementation (solving for velocity increment)
      timeint_gen_alpha,  // Peter's implementation (solving for acceleration increment)
      timeint_afgenalpha,
      timeint_bdf2,
    };

    /// initial field
    enum InitialField
    {
      initfield_zero_field,
      initfield_field_by_function,
      initfield_disturbed_field_from_function,
      initfield_flame_vortex_interaction,
      initfield_beltrami_flow,
      initfield_kim_moin_flow,
      initfield_bochev_test
    };

    /// initial field
    enum CalcError
    {
      no_error_calculation,
      beltrami_flow,
      channel2D,
      gravitation,
      shear_flow
    };

    /// initial field
    enum MeshTying
    {
      no_meshtying,
      condensed_smat,
      condensed_bmat,
      sps_coupled,
      sps_pc
    };

    //! physical turbulence models
    enum TurbModelAction
    {
      no_model,
      smagorinsky,
      smagorinsky_with_van_Driest_damping,
      dynamic_smagorinsky,
      scale_similarity,
      mixed_scale_similarity_eddy_viscosity_model
    };

    //! options for fine-scale subgrid viscosity
    enum FineSubgridVisc
    {
      no_fssgv,
      smagorinsky_all,
      smagorinsky_small
    };

    //! linearisation actions recognized by fluid3 (genalpha implementation)
    enum LinearisationAction
    {
      no_linearisation,
      fixed_point_like,
      Newton          ,
      minimal
    };

    //---------------------------------------------
    // some enums which could also be defined here
    //---------------------------------------------
    // more time integration schemes
    // linearization strategy convective term
    // type of stabilization strategy (SUPG,PSPG,CSTAB,...)

  }  // namespace FLUID

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_COMBUST_H
#endif  // #ifdef CCADISCRET
