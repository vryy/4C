/*----------------------------------------------------------------------*/
/*!
\file inpar_ale.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_ALE_H
#define INPAR_ALE_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for ALE
  namespace ALE
  {

  /// possible types of moving boundary simulation
  enum AleDynamic
  {
    // old types
    none,
    classic_lin, // will die
    incr_lin, // will die
    laplace,
    springs,
    springs_fixed_ref, // will die
    // new type
    solid,
  };
  /// Handling of non-converged nonlinear solver
  enum DivContAct
  {
    divcont_stop,               ///< abort simulation
    divcont_continue           ///< continue nevertheless
  };


  }  // namespace ALE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ALE_H
