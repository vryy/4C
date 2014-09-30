/*!-------------------------------------------------------------------------*
\file xfluid_functions.H

\brief Managing and evaluating of spatial functions for xfluid cases

<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*----------------------------------------------------------------------*/

#include "../drt_lib/drt_function.H"

#ifndef XFLUID_FUNCTION_H
#define XFLUID_FUNCTION_H

namespace DRT
{
  class Discretization;

  namespace UTILS
  {

    /// special implementation level set cut utilizing xfluid
    class GerstenbergerForwardfacingStep : public Function
    {
    public:

      /// ctor
      GerstenbergerForwardfacingStep();

      /// evaluate function at given position in space
      double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    };

  }
}

#endif
