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

    /// type of stabilization parameter definition
    //TODO:
    enum TauType
    {
      tautype_franca_barrenechea_valentin_wall,   // combustion
      tautype_bazilevs,                           // combustion
      smoothed_franca_barrenechea_valentin_wall,  // Fluid3
      franca_barrenechea_valentin_wall         ,  // Fluid3
      franca_barrenechea_valentin_codina       ,  // Fluid3
      fbvw_gradient_based_hk                   ,  // Fluid3
      codina                                   ,  // Fluid3
      bazilevs                                 ,  // Fluid3
      fbvw_wo_dt                               ,  // Fluid3
      tau_not_defined
    };

    /// time integration scheme
    enum TimeIntegrationScheme
    {
      timeint_stationary,
      timeint_one_step_theta
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
