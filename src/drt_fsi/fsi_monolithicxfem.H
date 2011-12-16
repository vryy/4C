#ifdef CCADISCRET

#ifndef FSI_MONOLITHICXFEM_H
#define FSI_MONOLITHICXFEM_H

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_fluid_xfem.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_io/io_ostream0.H"

#include "fsi_overlapprec.H"
#include "fsi_matrixtransform.H"
#include "fsi_statustest.H"
#include "fsi_debugwriter.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

#include "fsi_monolithicinterface.H"

namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;
}

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
    explicit MonolithicBaseXFEM(const Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~MonolithicBaseXFEM();

    /// read restart data
    void ReadRestart(int step);

    /// fluid field solver
    ADAPTER::FluidXFEM& FluidField() { return fluidfield_; }
    /// return instance of fluid time integration class
    const ADAPTER::FluidXFEM& FluidField() const { return fluidfield_; }

    const ADAPTER::Coupling& StructureFluidCoupling() const { return coupsf_; }



  protected:

    //! @name Time loop building blocks

    /// start a new time step
    virtual void PrepareTimeStep();

    /// take current results for converged and save for next time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput()
    {
      StructureField().PrepareOutput();
    }

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

    //!< a singular stream to stdout in parallel environment
    IO::Ostream0 cout0_;

  };


  /// base class of XFEM monolithic FSI algorithms
  class MonolithicXFEM : public MonolithicBaseXFEM
  {
  public:
    explicit MonolithicXFEM(const Epetra_Comm& comm);

    /// outer level FSI time loop
    void Timeloop();

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    /// Extract initial guess from fields
    virtual void InitialGuess();

    /// extractor to communicate between full monolithic map and block maps
    const MonolithicXFEMExtractor& Extractor() const {return extractor_;};

    /// full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const { return extractor_.FullMap(); }

    /// combined DBC map
    Teuchos::RCP<Epetra_Map> CombinedDBCMap();

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

    bool Converged();

    void LinearSolve();

    void VectorRescue(
        Teuchos::RCP<Epetra_Vector>& oldv
        ) const;

    MonolithicXFEMExtractor extractor_;

    Teuchos::RCP<Epetra_Vector> StructToFluid(const Teuchos::RCP<const Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> FluidToStruct(const Teuchos::RCP<const Epetra_Vector> iv) const;

    /// residual
    Teuchos::RCP<Epetra_Vector> rhs_;
    /// increment between timestep n and n+1
    Teuchos::RCP<Epetra_Vector> stepinc_;
    /// increment between Newton steps i and i+1
    Teuchos::RCP<Epetra_Vector> iterinc_;
    /// system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;
  };

}

#endif
#endif
