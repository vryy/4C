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
#ifdef CCADISCRET
#ifndef INPAR_SCATRA_H
#define INPAR_SCATRA_H

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
      velocity_Navier_Stokes
    };

    /// initial field for scalar transport problem
    enum InitialField
    {
      initfield_zero_field,
      initfield_field_by_function,
      initfield_field_by_condition,
      initfield_disturbed_field_by_function,
      initfield_DISCONTPV_1D,
      initfield_FVI_FERECHPRO,
      initfield_RAYTAYMIXFRAC
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
      calcerror_Kwok_Wu
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
      tau_franca_valentin,
      tau_bazilevs,
      tau_exact_1d,
      tau_zero
    };

    /// possible options for all-scale subgrid diffusivity
    enum KartType
    {
      kart_artificial,
      kart_hughes,
      kart_tezduyar,
      kart_docarmo,
      kart_almeida
    };

    /// possible options for fine-scale subgrid diffusivity
    enum FSSUGRDIFF
    {
      fssugrdiff_no,
      fssugrdiff_artificial,
      fssugrdiff_transfer_from_fluid
    };

    /// parameters for flux calculation of scalar transport problems
    enum FluxType
    {
      flux_no,
      flux_total_domain,
      flux_diffusive_domain,
      flux_total_boundary,
      flux_diffusive_boundary
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

    /// type of scalar transport problem
    enum ScaTraType
    {
      scatratype_standard,
      scatratype_levelset
    };

    /// type of reinitialization strategy for level set function
    enum ReinitializationAction
    {
      reinitaction_none,
      reinitaction_directdistance,
      reinitaction_sussman,
      reinitaction_interfaceprojection,
      reinitaction_function,
      reinitaction_signeddistancefunction
    };

    /// type of mass calculation strategy for level set function
    enum MassCalculation
    {
      masscalc_none,
      masscalc_squares,
      masscalc_interpolated
    };

  }  // namespace SCATRA

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_SCATRA_H
#endif  // #ifdef CCADISCRET
