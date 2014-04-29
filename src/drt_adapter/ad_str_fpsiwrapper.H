/*----------------------------------------------------------------------*/
/*!
\file ad_str_fpsiwrapper.H

\brief Structural adapter for FPSI problems containing the interface
       and methods dependent on the interface

<pre>
Maintainer: Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 - 15249
</pre>
*/

#ifndef FPSI_STRUCTURE_WRAPPER_H
#define FPSI_STRUCTURE_WRAPPER_H

#include "ad_str_fsiwrapper.H"

namespace ADAPTER
{

  class FPSIStructureWrapper : public FSIStructureWrapper
  {
  public:
    /// constructor
    explicit FPSIStructureWrapper(Teuchos::RCP<Structure> structure);

    /// Destructor
    virtual ~FPSIStructureWrapper() { ; }

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn(bool FPSI);

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp(bool FPSI);
  };
}
#endif
