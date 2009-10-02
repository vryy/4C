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

    /// type of calar transport problem
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
