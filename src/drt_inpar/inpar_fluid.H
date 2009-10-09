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
    // remark: more definitions to be added in the future
    enum TauType
    {
      tautype_franca_barrenechea_valentin_wall,
      tautype_bazilevs
    };
    //---------------------------------------------
    // some enums which could also be defined here
    //---------------------------------------------
    // time integration scheme
    // linearization strategy convective term
    // type of stabilization strategy (SUPG,PSPG,CSTAB,...)

  }  // namespace FLUID

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_COMBUST_H
#endif  // #ifdef CCADISCRET
