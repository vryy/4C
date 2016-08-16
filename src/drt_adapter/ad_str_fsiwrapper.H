/*----------------------------------------------------------------------*/
/*!
\file ad_str_fsiwrapper.H

\brief Structural adapter for FSI problems containing the interface
       and methods dependent on the interface

\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237

\level 1
*/

#ifndef FSI_STRUCTURE_WRAPPER_H
#define FSI_STRUCTURE_WRAPPER_H

#include "ad_str_wrapper.H"
#include "../drt_fsi/fsi_str_model_evaluator_partitioned.H"

// forward declarations
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}

namespace STR
{
  namespace MODELEVALUATOR
  {
    class PartitionedFSI;
  }
}


namespace ADAPTER
{

  class FSIStructureWrapper : public StructureWrapper
  {
  public:

    /// constructor
    explicit FSIStructureWrapper(Teuchos::RCP<Structure> structure);

    /// Destructor
    virtual ~FSIStructureWrapper() { ; }

    /// communication object at the interface
    virtual Teuchos::RCP<const STR::AUX::MapExtractor> Interface() const
    {
      return interface_;
    }

    /// switch structure field to block matrix in fsi simulations
    virtual void UseBlockMatrix();

    /// linear structure solve with just a interface load
    ///
    /// The very special solve done in steepest descent relaxation
    /// calculation (and matrix free Newton Krylov).
    ///
    /// \note Can only be called after a valid structural solve.
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce);

    /// @name Extract interface values

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn();

    /// extract interface displacements at \f$t_{n+1}\f$
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

    /// remove as soon as new structure is fully usable ! todo
    /// only 4 nightly tests use this method:
    /// fsi_dc3D_part_ait_ga_ost_xwall (no solidsh8 possible yet)
    /// fsi_ow3D_mtr_drt (no solidsh8 possible yet)
    /// crfsi_inclDomain_suppFluid_interfaceBuild (crack specific model evaluator is missing)
    /// constr2D_fsi (newtonlinuzawa not implemented; but really needed ?)
    void ApplyInterfaceForcesTemporaryDeprecated(Teuchos::RCP<Epetra_Vector> iforce);

    /// rebuild FSI interface from structure side
    virtual void RebuildInterface();

    /// set pointer to model evaluator
    void SetModelEvaluatorPtr(Teuchos::RCP<STR::MODELEVALUATOR::PartitionedFSI> me)
    {
      fsi_model_evaluator_=me;
      return;
    }

  protected:

    /// the interface map setup for interface <-> full translation
    Teuchos::RCP<STR::AUX::MapExtractor> interface_;

    /// predictor type
    int predictor_;

    /// The structural model evaluator object.
    /// Your FSI algorithm calls methods in this adapter.
    /// If this method is related to the structural field,
    /// a corresponding method in the model evaluator may be
    /// called, if necessary.
    /// See e.g. \ref ADAPTER::FSIStructureWrapper::RelaxationSolve()
    Teuchos::RCP<STR::MODELEVALUATOR::PartitionedFSI> fsi_model_evaluator_;

    /// access the fsi model evaluator
    Teuchos::RCP<STR::MODELEVALUATOR::PartitionedFSI> FSIModelEvaluator(){return fsi_model_evaluator_;};

  }; // class FSIStructureWrapper
} // namespace ADAPTER

#endif
