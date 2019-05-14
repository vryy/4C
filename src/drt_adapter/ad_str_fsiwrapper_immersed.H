/*----------------------------------------------------------------------*/
/*!

\brief Structural adapter for Immersed and Immersed + ALE FSI problems containing the interface
       and methods dependent on the interface

\level 2

\maintainer Jonas Eichinger

*/

#ifndef FSI_STRUCTURE_WRAPPER_IMMERSED_H
#define FSI_STRUCTURE_WRAPPER_IMMERSED_H

#include "ad_str_fpsiwrapper.H"
#include "../drt_structure_new/str_dbc.H"

// forward declarations
namespace STR
{
  class Dbc;

  namespace AUX
  {
    class MapExtractor;
  }
}  // namespace STR


namespace ADAPTER
{
  class FSIStructureWrapperImmersed : public FPSIStructureWrapper
  {
   public:
    /// constructor
    explicit FSIStructureWrapperImmersed(Teuchos::RCP<Structure> structure);

    /// destructor
    virtual ~FSIStructureWrapperImmersed() { ; }

    /// return the combined interface (matching fsi + immersed)
    const Teuchos::RCP<LINALG::MapExtractor> CombinedInterface() { return combinedinterface_; };

    /// extract interface displacements at \f$t_{n+1}\f$ of immersed interface
    virtual Teuchos::RCP<Epetra_Vector> ExtractImmersedInterfaceDispnp();

    /// extract interface displacements at \f$t_{n+1}\f$ of immersed interface + fsi interface
    virtual Teuchos::RCP<Epetra_Vector> ExtractFullInterfaceDispnp();

    /// Predictor for interface displacements at immersed interface
    Teuchos::RCP<Epetra_Vector> PredictImmersedInterfaceDispnp();

    /// Predictor for interface displacements at immersed and fsi interface
    Teuchos::RCP<Epetra_Vector> PredictFullInterfaceDispnp();

    /// Get mutable reference to DBC object
    STR::Dbc& GetMutableDBC();

    /// expand dirichlet bc map
    /// old struct. time integration version
    virtual void AddDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoadd);

    /// contract dirichlet bc map
    /// old struct. time integration version
    virtual void RemoveDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoremove);

    /// set the state of the nox group and the global state data container
    void SetState(const Teuchos::RCP<Epetra_Vector>& x);

    /// @name Apply interface forces

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    void ApplyImmersedInterfaceForces(
        Teuchos::RCP<Epetra_Vector> iforce_fsi, Teuchos::RCP<Epetra_Vector> iforce_immersed);

    /*!
      \brief Write extra output for specified step and time.
             Useful if you want to write output every iteration in partitioned schemes.
             If no step and time is provided, standard Output of structure field is invoked.

      \param forced_writerestart (in) : Force to write restart
      \param step (in) : Pseudo-step for which extra output is written
      \param time (in) : Pseudo-time for which extra output is written

      \note This is a pure DEBUG functionality. Originally used in immersed method development.

      \warning This method partly re-implements redundantly few lines of the common structure
      Output() routine. \return void
    */
    virtual void Output(
        bool forced_writerestart = false, const int step = -1, const double time = -1.0);

   protected:
    /// the interface map setup for immersed interface <-> fsi interface distinction
    Teuchos::RCP<LINALG::MapExtractor> combinedinterface_;

    /// combined matching FSI - IMMERSED interface map
    Teuchos::RCP<Epetra_Map> combinedmap_;
  };
}  // namespace ADAPTER

#endif
