/*----------------------------------------------------------------------*/
/*!
\file inpar_ale.H

<pre>
 Maintainer: Matthias Mayr
             mayr@mhpc.mw.tum.de
             089 - 289-10362
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
    none, ///< default
    laplace, ///< Laplacian smoothing
    springs, ///< use a spring analogy
    solid, ///< pseudo-structure approach
  };
  /// Handling of non-converged nonlinear solver
  enum DivContAct
  {
    divcont_stop, ///< abort simulation
    divcont_continue ///< continue nevertheless
  };


  }  // namespace ALE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ALE_H
