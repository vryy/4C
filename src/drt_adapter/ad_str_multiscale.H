/*----------------------------------------------------------------------*/
/*!
\file ad_str_multiscale.H

\brief Wrapper for the structural time loop for multiscale problems

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_MULTISCALE_H
#define ADAPTER_STRUCTURE_MULTISCALE_H

#include "ad_str_timeloop.H"

namespace ADAPTER
{

  /// wrap for multi scale specific implementation
  class StructureMultiScale : public StructureTimeLoop
  {
  public:

    /// constructor
    explicit StructureMultiScale(Teuchos::RCP<Structure> structure) : StructureTimeLoop(structure) {}

    /// stop possible supporting procs after the time loop
    virtual void PostTimeLoop();

  };

}

#endif
