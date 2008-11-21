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

  }  // namespace SCATRA

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_SCATRA_H
#endif  // #ifdef CCADISCRET
