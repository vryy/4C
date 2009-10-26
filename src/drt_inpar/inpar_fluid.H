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
    /// type of stabilization parameter definition
    enum TauType
    {
      tautype_franca_barrenechea_valentin_wall,
      tautype_bazilevs
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
