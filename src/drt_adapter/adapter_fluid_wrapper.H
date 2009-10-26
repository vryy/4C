#ifndef ADAPTER_FLUID_WRAPPER_H
#define ADAPTER_FLUID_WRAPPER_H

#include "adapter_fluid.H"

namespace ADAPTER
{

  /// Just wrap, do nothing new, meant to be derived from
  class FluidWrapper : public Fluid
  {
  public:

    explicit FluidWrapper(Teuchos::RCP<Fluid> fluid) : fluid_(fluid) {}

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess()
    { return fluid_->InitialGuess(); }
    virtual Teuchos::RCP<const Epetra_Vector> RHS()
    { return fluid_->RHS(); }
    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual()
    { return fluid_->TrueResidual(); }
    virtual Teuchos::RCP<const Epetra_Vector> Velnp()
    { return fluid_->Velnp(); }
    virtual Teuchos::RCP<const Epetra_Vector> Velaf()
    { return fluid_->Velaf(); }
    virtual Teuchos::RCP<const Epetra_Vector> Veln()
    { return fluid_->Veln(); }
    virtual Teuchos::RCP<const Epetra_Vector> Accam()
    { return fluid_->Accam(); }
    virtual Teuchos::RCP<const Epetra_Vector> Hist()
    { return fluid_->Hist(); }
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp()
    { return fluid_->Dispnp(); }
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel()
    { return fluid_->ConvectiveVel(); }
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap()
    { return fluid_->DofRowMap(); }
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
    { return fluid_->SystemMatrix(); }
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    { return fluid_->BlockSystemMatrix(); }
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives()
    { return fluid_->ShapeDerivatives(); }
    virtual Teuchos::RCP<DRT::Discretization> Discretization()
    { return fluid_->Discretization(); }
    virtual Teuchos::RCP<const DRT::DofSet> DofSet()
    { return fluid_->DofSet(); }
    virtual void SetInitialFlowField(int whichinitialfield,int startfuncno)
    { return fluid_->SetInitialFlowField(whichinitialfield, startfuncno); }
    virtual void SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,
                                   RCP<const Epetra_Vector> scalaram,
                                   RCP<const Epetra_Vector> scalardtam,
                                   const double             thermpressaf,
                                   const double             thermpressam,
                                   const double             thermpressdtam,
                                   const int                numscal)
    { return fluid_->SetIterLomaFields(scalaraf, scalaram, scalardtam, thermpressaf,
                                       thermpressam, thermpressdtam, numscal); }
    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,
                                   const double             thermpressnp,
                                   RCP<const Epetra_Vector> scatraresidual,
                                   const int                numscal)
    { return fluid_->SetTimeLomaFields(scalarnp, thermpressnp, scatraresidual, numscal); }
    virtual void PrepareTimeStep()
    { return fluid_->PrepareTimeStep(); }
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
    { return fluid_->Evaluate(); }
    virtual void Update()
    { return fluid_->Update(); }
    virtual void StatisticsAndOutput()
    { return fluid_->StatisticsAndOutput(); }
    virtual void Output()
    { return fluid_->Output(); }
    virtual void ReadRestart(int step)
    { return fluid_->ReadRestart(); }
    virtual double Time() const
    { return fluid_->Time(); }
    virtual int Step() const
    { return fluid_->Step(); }
    virtual void LiftDrag()
    { return fluid_->LiftDrag(); }
    virtual void NonlinearSolve()
    { return fluid_->NonlinearSolve(); }
    virtual void Predictor()
    { return fluid_->Predictor(); }
    virtual void MultiCorrector()
    { return fluid_->MultiCorrector(); }
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
    { return fluid_->RelaxationSolve(ivel); }
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver()
    { return fluid_->LinearSolver(); }
    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap()
    { return fluid_->InnerVelocityRowMap(); }
    virtual Teuchos::RCP<const Epetra_Map>    VelocityRowMap()
    { return fluid_->VelocityRowMap(); }
    virtual Teuchos::RCP<const Epetra_Map> PressureRowMap()
    { return fluid_->PressureRowMap(); }
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
    { return fluid_->SetMeshMap(mm); }
    virtual double ResidualScaling() const
    { return fluid_->ResidualScaling(); }
    virtual double TimeScaling() const
    { return fluid_->TimeScaling(); }
    virtual const FLD::UTILS::MapExtractor& Interface() const
    { return fluid_->Interface(); }
    virtual const FLUID_TIMEINTTYPE TimIntScheme() const
    { return fluid_->TimIntScheme(); }
    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
    { return fluid_->ExtractVelocityPart(velpres); }
    virtual void ImportDiscretization(Teuchos::RCP<const DRT::Discretization> importdis)
    {
      dserror("This function is not implemented for this type of fluid! Only available for combustion!");
    };
    virtual void ImportInterface(const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& importinterface)
    {
      dserror("This function is not implemented for this type of fluid! Only available for combustion!");
    };
    virtual void ImportFlameFront(Teuchos::RCP<const COMBUST::FlameFront> importflamefront)
    {
      dserror("This function is not implemented for this type of fluid! Only available for combustion!");
    };
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
    { return fluid_->ApplyInterfaceVelocities(ivel); }
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel,Teuchos::RCP<Epetra_Vector> iforce)
    { return fluid_->ApplyInterfaceRobinValue(ivel,iforce); }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity()
    { return fluid_->ExtractInterfaceFluidVelocity(); }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln() = 0;
    { return fluid_->ExtractInterfaceVeln(); }
    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln()
    { return fluid_->ExtractFreeSurfaceVeln(); }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces()
    { return fluid_->ExtractInterfaceForces(); }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin()
    { return fluid_->ExtractInterfaceForcesRobin(); }
    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp)
    { return fluid_->ApplyMeshDisplacement(fluiddisp); }
    virtual void ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel)
    { return fluid_->ApplyMeshVelocity(gridvel); }
    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
    { return fluid_->DisplacementToVelocity(fcx); }
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
    { return fluid_->VelocityToDisplacement(fcx); }
    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
    { return fluid_->FreeSurfDisplacementToVelocity(fcx); }
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
    { return fluid_->FreeSurfVelocityToDisplacement(fcx); }
    virtual int Itemax() const
    { return fluid_->Itemax(); }
    virtual void SetItemax(int itemax)
    { return fluid_->SetItemax(itemax); }
    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape()
    { return fluid_->IntegrateInterfaceShape(); }
    virtual void UseBlockMatrix(bool splitmatrix)
    { return fluid_->UseBlockMatrix(splitmatrix); }
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest()
    { return fluid_->CreateFieldTest(); }


  protected:

    Teuchos::RCP<Structure> fluid_;
  }
}

#endif
