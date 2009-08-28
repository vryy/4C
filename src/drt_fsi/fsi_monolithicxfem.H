#ifdef CCADISCRET

#ifndef FSI_MONOLITHICXFEM_H
#define FSI_MONOLITHICXFEM_H

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_fluid_xfem.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_lib/linalg_mapextractor.H"

#include "fsi_overlapprec.H"
#include "fsi_matrixtransform.H"
#include "fsi_statustest.H"
#include "fsi_debugwriter.H"
#include "fsi_monolithicinterface.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

namespace FSI
{

  class MonolithicXFEMExtractor : public LINALG::MultiMapExtractor
  {
  public:

    MAP_EXTRACTOR_VECTOR_METHODS(Structure,0)
    MAP_EXTRACTOR_VECTOR_METHODS(Fluid,    1)

  };

  /// base class for monolithic FSI with xfem
  class MonolithicBaseXFEM : public ADAPTER::AlgorithmBase,
                             public ADAPTER::StructureBaseAlgorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit MonolithicBaseXFEM(Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~MonolithicBaseXFEM();

    /// read restart data
    void ReadRestart(int step);

    /// fluid field solver
    ADAPTER::FluidXFEM& FluidField() { return fluidfield_; }
    /// return instance of fluid time integration class
    const ADAPTER::FluidXFEM& FluidField() const { return fluidfield_; }

    ADAPTER::Coupling& StructureFluidCoupling() { return coupsf_; }

  protected:

    //! @name Time loop building blocks

    /// start a new time step
    virtual void PrepareTimeStep();

    /// take current results for converged and save for next time step
    virtual void Update();

    /// write output
    virtual void Output();

  protected:

    /// coupling of structure and fluid at the interface
    ADAPTER::Coupling               coupsf_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks
    UTILS::MatrixRowColTransform sggtransform_;
    UTILS::MatrixRowTransform    sgitransform_;
    UTILS::MatrixColTransform    sigtransform_;
    ///@}

    /// the very xfem fluid
    ADAPTER::FluidXFEM              fluidfield_;
  };


  /// monolithic FSI with xfem
  class MonolithicXFEM : public MonolithicBaseXFEM
  {
  public:
    explicit MonolithicXFEM(Epetra_Comm& comm);

    /// outer level FSI time loop
    void Timeloop();

  private:

    void Newton();
    void Evaluate();
    void SetupRHS();
    bool ConverganceTest();
    void SetupSystemMatrix();
    void LinearSolve();

    Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv);
    Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv);

    Teuchos::RCP<MonolithicXFEMExtractor> extractor_;

    Teuchos::RCP<Epetra_Vector> dispincsum_;
    Teuchos::RCP<Epetra_Vector> rhs_;
    Teuchos::RCP<Epetra_Vector> x_;
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> matrix_;
  };

}

#endif
#endif
