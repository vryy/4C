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
    classic_lin,
    incr_lin,
    laplace,
    springs,
    springs_fixed_ref
  };

  }  // namespace ALE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ALE_H
