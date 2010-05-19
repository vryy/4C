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

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

#include "fsi_monolithicinterface.H"
namespace FSI
{

  /// monolithic FSI with XFEM: map extractor
  class MonolithicXFEMExtractor : public LINALG::MultiMapExtractor
  {
  public:

    MAP_EXTRACTOR_VECTOR_METHODS(StructureInterior,0)
    MAP_EXTRACTOR_VECTOR_METHODS(StructureBoundary,1)
    MAP_EXTRACTOR_VECTOR_METHODS(Fluid            ,2)

  };

  /// base class for monolithic FSI with XFEM
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

    //@}



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


  /// base class of XFEM monolithic FSI algorithms
  class MonolithicXFEM : public MonolithicBaseXFEM
  {
  public:
    explicit MonolithicXFEM(Epetra_Comm& comm);

    /// outer level FSI time loop
    void Timeloop();


    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    /// extract the two field vectors from a given composed vector
    /*!
      x is the sum of all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
     */
    void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                             Teuchos::RCP<const Epetra_Vector>& sx,
                             Teuchos::RCP<const Epetra_Vector>& fx);


    /// extractor to communicate between full monolithic map and block maps
    const MonolithicXFEMExtractor& Extractor(){return extractor_;};

    /// full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const { return extractor_.FullMap(); }

    /// set full monolithic dof row map
    /*!
      A subclass calls this method (from its constructor) and thereby
      defines the number of blocks, their maps and the block order. The block
      maps must be row maps by themselves and must not contain identical GIDs.
     */
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    void SetupExtractor();

  private:

//    void SetupVector(
//        Epetra_Vector &f,
//        Teuchos::RCP<const Epetra_Vector> sv,
//        Teuchos::RCP<const Epetra_Vector> fv,
//        double fluidscale);

    void Newton();
    void Evaluate();
    /// setup composed right hand side from field solvers
    void SetupRHS();

    /// setup composed system matrix from field solvers
    void SetupSystemMatrix();

    bool ConverganceTest();

    void LinearSolve();

    MonolithicXFEMExtractor extractor_;

    Teuchos::RCP<Epetra_Vector> StructToFluid(const Teuchos::RCP<const Epetra_Vector> iv);
    Teuchos::RCP<Epetra_Vector> FluidToStruct(const Teuchos::RCP<const Epetra_Vector> iv);

    /// residual
    Teuchos::RCP<Epetra_Vector> rhs_;
    /// increment between timestep n and n+1
    Teuchos::RCP<Epetra_Vector> stepinc_;
    /// system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;
  };

}

#endif
#endif
