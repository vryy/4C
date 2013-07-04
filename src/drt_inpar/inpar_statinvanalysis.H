/*----------------------------------------------------------------------*/
/*!
\file inpar_statinvanalysis.H

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_STATINVANALYSIS_H
#define INPAR_STATINVANALYSIS_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for Inverse Analysis
  namespace STR
  {
    enum StatInvMatParametrization
    {
      stat_inv_mp_none,                  // nothing defined
      stat_inv_mp_smoothkernel,          // smooth kernel approximation
      stat_inv_mp_elementwise            // parameters per element
    };

    /// Type of Inverse Analysis
    enum StatInvAnalysisType
    {
      stat_inv_none,          // no
      stat_inv_mc,            // statistical inverse analysis using montecarlo
      stat_inv_graddesc       // deterministic inverse analysis using gradient descent
    };

  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATINVANALYSIS_H
