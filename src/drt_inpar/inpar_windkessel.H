/*----------------------------------------------------------------------*/
/*!
\file inpar_windkessel.H

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef INPAR_WINDKESSEL_H_
#define INPAR_WINDKESSEL_H_

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

namespace INPAR
{

namespace WINDKESSEL
{

  /// possible windkessel-structural solvers
  enum WindkSolveAlgo
  {
    windksolve_direct, ///< build monolithic windkessel-structural system
    windksolve_simple  ///< use simple preconditioner for iterative solve
  };

  /// set the windkessel parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  /// set specific windkessel conditions
  void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

}// namespace WINDKESSEL
}// namespace INPAR
/*----------------------------------------------------------------------*/
#endif /* INPAR_WINDKESSEL_H_ */
