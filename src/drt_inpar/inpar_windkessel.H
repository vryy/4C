/*----------------------------------------------------------------------*/
/*!
\file inpar_windkesseö.H

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

}// namespace WINDKESSEL
}// namespace INPAR
/*----------------------------------------------------------------------*/
#endif /* INPAR_WINDKESSEL_H_ */
