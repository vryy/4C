/*----------------------------------------------------------------------*/
/*!
\file ad_str_statmech.H

\brief Wrapper for the time loop of statmech problems

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_STATMECH_H
#define ADAPTER_STRUCTURE_STATMECH_H

#include "ad_str_timeloop.H"

namespace ADAPTER
{

  /// wrap for statmech specific implementation
  class StructureStatMech : public StructureTimeLoop
  {
  public:

    /// constructor
    explicit StructureStatMech(Teuchos::RCP<Structure> structure) : StructureTimeLoop(structure) {}

    /// statmech specific time loop
    virtual int Integrate();

  };

}

#endif
