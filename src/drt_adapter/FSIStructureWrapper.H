#ifndef FSI_STRUCTURE_WRAPPER_H
#define FSI_STRUCTURE_WRAPPER_H

#include "adapter_structure_wrapper.H"


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

  /// Calculate increments from absolute values
  class FSIStructureWrapper : public StructureWrapper
  {
  public:

    explicit FSIStructureWrapper(Teuchos::RCP<Structure> structure);


    /// communication object at the interface
    virtual const STR::AUX::MapExtractor& Interface() const
    {
     return *interface_;
    }

    //! switch structure field to block matrix in fsi simulations
    virtual void UseBlockMatrix();

    /// linear structure solve with just a interface load
    ///
    /// The very special solve done in steepest descent relaxation
    /// calculation (and matrix free Newton Krylov).
    ///
    /// \note Can only be called after a valid structural solve.
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce);

    /// @name Extract interface values

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn();

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp();

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp();

    /// @name Apply interface forces

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce);


  private:

    /// the interface map setup for interface <-> full translation
    Teuchos::RCP<STR::AUX::MapExtractor> interface_;

    /// predictor type
    int predictor_;

  };
}

#endif
