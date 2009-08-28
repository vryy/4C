#ifndef ADAPTER_STRUCTURE_WRAPPER_H
#define ADAPTER_STRUCTURE_WRAPPER_H

#include "adapter_structure.H"

namespace ADAPTER
{

  /// Just wrap, do nothing new, meant to be derived from
  class StructureWrapper : public Structure
  {
  public:

    explicit StructureWrapper(Teuchos::RCP<Structure> structure) : structure_(structure) {}

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess()
    { return structure_->InitialGuess(); }

    virtual Teuchos::RCP<const Epetra_Vector> RHS()
    { return structure_->RHS(); }
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp()
    { return structure_->Dispnp(); }
    virtual Teuchos::RCP<const Epetra_Vector> Dispn()
    { return structure_->Dispn(); }

    virtual Teuchos::RCP<const Epetra_Map> DofRowMap()
    { return structure_->DofRowMap(); }
    virtual const Epetra_Map& DomainMap()
    { return structure_->DomainMap(); }
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
    { return structure_->SystemMatrix(); }
    virtual Teuchos::RCP<DRT::Discretization> Discretization()
    { return structure_->Discretization(); }
    virtual bool HaveConstraint()
    { return structure_->HaveConstraint(); }
    virtual const RCP<const UTILS::ConstrManager> GetConstraintManager()
    { return structure_->GetConstraintManager(); }
    virtual const RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    { return structure_->GetDBCMapExtractor(); }

    virtual const LINALG::MapExtractor& Interface() const
    { return structure_->Interface(); }
    virtual double GetTime() const
    { return structure_->GetTime(); }
    virtual double GetTimeEnd() const
    { return structure_->GetTimeEnd(); }
    virtual double GetTimeStepSize() const
    { return structure_->GetTimeStepSize(); }
    virtual int GetTimeStep() const
    { return structure_->GetTimeStep(); }
    virtual int GetTimeNumStep() const
    { return structure_->GetTimeNumStep(); }

    virtual void PrepareTimeStep()
    { structure_->PrepareTimeStep(); }
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
    { structure_->Evaluate(disp); }
    virtual void Update()
    { structure_->Update(); }
    virtual void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr)
    { structure_->UpdateIterIncrConstr(lagrincr); }


    virtual void Output()
    { return structure_->Output(); }
    virtual void ReadRestart(int step)
    { return structure_->ReadRestart(step); }

    virtual void Solve()
    { return structure_->Solve(); }
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce)
    { return structure_->RelaxationSolve(iforce); }
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver()
    { return structure_->LinearSolver(); }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn()
    { return structure_->ExtractInterfaceDispn(); }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp()
    { return structure_->ExtractInterfaceDispnp(); }
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp()
    { return structure_->PredictInterfaceDispnp(); }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces()
    { return structure_->ExtractInterfaceForces(); }

    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce)
    { structure_->ApplyInterfaceForces(iforce); }
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> iforce,
                                          Teuchos::RCP<Epetra_Vector> ifluidvel)
    {  structure_->ApplyInterfaceRobinValue(iforce,ifluidvel); }
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest()
    { return structure_->CreateFieldTest(); }

  protected:

    Teuchos::RCP<Structure> structure_;
  };


  /// Calculate increments from absolute values
  class StructureNOXCorrectionWrapper : public StructureWrapper
  {
  public:

    explicit StructureNOXCorrectionWrapper(Teuchos::RCP<Structure> structure)
      : StructureWrapper(structure) {}

    virtual void PrepareTimeStep();
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disp);

  private:

    /// sum of displacement vectors already applied,
    /// i.e. the incremental displacements
    Teuchos::RCP<Epetra_Vector> disinc_;
  };
}

#endif
