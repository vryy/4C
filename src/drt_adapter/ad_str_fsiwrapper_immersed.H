/*----------------------------------------------------------------------*/
/*!
\file ad_str_fsiwrapper_immersed.H

\brief Structural adapter for Immersed and Immersed + ALE FSI problems containing the interface
       and methods dependent on the interface

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

#ifndef FSI_STRUCTURE_WRAPPER_IMMERSED_H
#define FSI_STRUCTURE_WRAPPER_IMMERSED_H

#include "ad_str_fpsiwrapper.H"

// forward declarations
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}


namespace ADAPTER
{

  class FSIStructureWrapperImmersed : public FPSIStructureWrapper
  {
  public:

    /// constructor
    explicit FSIStructureWrapperImmersed(Teuchos::RCP<Structure> structure);

    /// Destructor
    virtual ~FSIStructureWrapperImmersed() { ; }

    const Teuchos::RCP<LINALG::MapExtractor> CombinedInterface(){return combinedinterface_;};

    /// extract interface displacements at \f$t_{n+1}\f$ of immersed interface
    virtual Teuchos::RCP<Epetra_Vector> ExtractImmersedInterfaceDispnp();

    /// extract interface displacements at \f$t_{n+1}\f$ of immersed interface + fsi interface
    virtual Teuchos::RCP<Epetra_Vector> ExtractFullInterfaceDispnp();

    /// Predictor for interface displacements at immersed interface
    Teuchos::RCP<Epetra_Vector> PredictImmersedInterfaceDispnp();

    /// Predictor for interface displacements at immersed and fsi interface
    Teuchos::RCP<Epetra_Vector> PredictFullInterfaceDispnp();

    /// @name Apply interface forces

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    /// special routine needed in immersed+ale fsi
    void ApplyImmersedInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce_fsi,
                                      Teuchos::RCP<Epetra_Vector> iforce_immersed);

    /// output results
    virtual void Output(bool forced_writerestart = false, const int step = -1, const double time = -1.0);

  protected:

    /// the interface map setup for immersed interface <-> fsi interface distinction
    Teuchos::RCP<LINALG::MapExtractor> combinedinterface_;

    /// combined FSI - IMMERSED interface map
    Teuchos::RCP<Epetra_Map> combinedmap_;

  };
}

#endif
