/*----------------------------------------------------------------------*/
/*!
\file inpar_scatra.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_SCATRA_H
#define INPAR_SCATRA_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for SCATRA
  namespace SCATRA
  {

    /// time integration schemes for scalar transport problems
    enum TimeIntegrationScheme
    {
      timeint_stationary,
      timeint_one_step_theta,
      timeint_bdf2,
      timeint_gen_alpha
    };

    /// type of solution procedures for scalar transport problems
    enum SolverType
    {
      solvertype_linear_full,
      solvertype_linear_incremental,
      solvertype_nonlinear
    };

    /// type of convective velocity field
    enum VelocityField
    {
      velocity_zero,
      velocity_function,
      velocity_function_and_curve,
      velocity_Navier_Stokes
    };

    /// initial field for scalar transport problem
    enum InitialField
    {
      initfield_zero_field,
      initfield_field_by_function,
      initfield_field_by_condition,
      initfield_disturbed_field_by_function,
      initfield_discontprogvar_1D,
      initfield_flame_vortex_interaction,
      initfield_raytaymixfrac,
      initfield_Lshapeddomain,
      initfield_facing_flame_fronts,
      initfield_oracles_flame,
      initialfield_forced_hit_high_Sc,
      initialfield_forced_hit_low_Sc
    };

    /// form of convective term
    enum ConvForm
    {
      convform_convective,
      convform_conservative
    };

    /// compute error compared to analytical solution
    enum CalcError
    {
      calcerror_no,
      calcerror_Kwok_Wu,
      calcerror_cylinder,
      calcerror_electroneutrality
    };

    /// possible types of stabilization
    enum StabType
    {
      stabtype_no_stabilization,
      stabtype_SUPG,
      stabtype_GLS,
      stabtype_USFEM
    };

    /// possible options for the stabilization parameter
    enum TauType
    {
      tau_taylor_hughes_zarins,
      tau_taylor_hughes_zarins_wo_dt,
      tau_franca_valentin,
      tau_franca_valentin_wo_dt,
      tau_shakib_hughes_codina,
      tau_shakib_hughes_codina_wo_dt,
      tau_codina,
      tau_codina_wo_dt,
      tau_franca_madureira_valentin,
      tau_franca_madureira_valentin_wo_dt,
      tau_exact_1d,
      tau_zero
    };

    /// possible options for characteristic element length
    enum CharEleLength
    {
      streamlength,
      volume_equivalent_diameter,
      root_of_volume
    };

    /// possible options for all-scale subgrid diffusivity
    enum AssgdType
    {
      assgd_artificial,
      assgd_hughes,
      assgd_tezduyar,
      assgd_tezduyar_wo_phizero,
      assgd_docarmo,
      assgd_almeida,
      assgd_lin_reinit,
      assgd_yzbeta,
      assgd_codina,
    };

    /// possible options for fine-scale subgrid diffusivity
    enum FSSUGRDIFF
    {
      fssugrdiff_no,
      fssugrdiff_artificial,
      fssugrdiff_smagorinsky_all,
      fssugrdiff_smagorinsky_small
    };

    /// parameters for flux calculation of scalar transport problems
    enum FluxType
    {
      flux_no,
      flux_total_domain,
      flux_diffusive_domain,
      flux_total_boundary,
      flux_diffusive_boundary,
      flux_convective_boundary
    };

    // this parameter selects the location where tau is evaluated
    enum EvalTau
    {
      evaltau_element_center,
      evaltau_integration_point
    };

    // this parameter selects the location where the material is evaluated
    enum EvalMat
    {
      evalmat_element_center,
      evalmat_integration_point
    };

    enum ScaTraType
    {
      scatratype_undefined,
      scatratype_condif,
      scatratype_loma,
      scatratype_elch,
      scatratype_levelset,
      scatratype_poro,
      scatratype_advreac,
      scatratype_pororeac,
      scatratype_cardio_monodomain
    };

    //! this parameter is closely connected to ScaTraType
    enum ImplType
    {
      impltype_undefined,
      impltype_std,
      impltype_loma,
      impltype_elch_NP,
      impltype_elch_diffcond,
      impltype_lsreinit,
      impltype_levelset,
      impltype_poro,
      impltype_advreac,
      impltype_pororeac,
      impltype_std_meshfree
    };

    /// type of method for improving consistency of stabilized methods
    enum Consistency
    {
      consistency_no,
      consistency_l2_projection_lumped
    };

    /// type of method for improving consistency of stabilized methods
    enum Kinetics
    {
      butler_volmer,
      butler_volmer_yang1997,
      tafel,
      linear,
      butler_volmer_newman,
      butler_volmer_bard,
      nernst,
      zero
    };

    /// type of method for improving consistency of stabilized methods
    enum Fluxeval
    {
      fluxeval_standard,
      fluxeval_alldof
    };

    enum Initpotcalc
    {
      initpotcalc_no,
      initpotcalc_yes
    };
  }  // namespace SCATRA

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_SCATRA_H
