/*----------------------------------------------------------------------*/
/*!
\file inpar_cavitation.H

\brief Input parameters for cavitation

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/
#ifndef INPAR_CAVITATION_H
#define INPAR_CAVITATION_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/

namespace INPAR
{

  /// Input parameters for cavitation
  namespace CAVITATION
  {

    /// Type of coupling strategy for cavitation problems
    enum CouplingStrategyOverFields
    {
      OneWay,
      TwoWayMomentum,
      TwoWayFull
    };


  }  // namespace CAVITATION

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CAVITATION_H
