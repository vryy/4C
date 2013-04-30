/*----------------------------------------------------------------------*/
/*!
\file inpar_crack.H

\brief Input parameters for crack propagation

<pre>
Maintainer: sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_CRACK_H
#define INPAR_CRACK_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for Inverse Analysis
  namespace CRACK
  {

    /// Type of crack modeling
    enum CrackModel
    {
      crack_none,          // no carack propagation
      crack_dczm,          // Discrete cohesive zone modeling
      crack_ddzm           // Discrete damage zone modeling
    };

    enum tractionSeparation
    {
      linear,
      trapezoidal,
      exponential,
      sinusoidal,
      ppr
    };

  }  // namespace CRACK

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif
