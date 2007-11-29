/*----------------------------------------------------------------------*/
/*!
\file drt_validparameters.H

\brief Setup of the list of valid input parameters

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef DRT_VALIDPARAMETERS_H
#define DRT_VALIDPARAMETERS_H

#ifdef __cplusplus /* stuff which is c++ specific */

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace DRT {

  /// construct list with all parameters and documentation
  Teuchos::RCP<const Teuchos::ParameterList> ValidParameters();

  /// set parameters read by our solvers
  void SetValidSolverParameters(Teuchos::ParameterList& list);

  /// set parameters for time adaptivity
  void SetValidTimeAdaptivityParameters(Teuchos::ParameterList& list);
}

#endif /* end of #ifdef __cplusplus */



/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C"
{
#endif

  /// print list of valid parameters with documentation
  void PrintValidParameters();

#ifdef __cplusplus
}
#endif

#endif
#endif
