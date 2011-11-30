#ifdef CCADISCRET

#ifndef FSI_FLUIDFLUIDMONOLITHICSTRUCTURESPLIT_H
#define FSI_FLUIDFLUIDMONOLITHICSTRUCTURESPLIT_H

#include "../drt_inpar/inpar_fsi.H"
#include "fsi_monolithic_nox.H"
#include "fsi_monolithic.H"
#include "fsi_overlapprec.H"
#include "fsi_matrixtransform.H"

namespace FSI
{

  /// monolithic Fluid-Fluid FSI algorithm Structure Split
  /*!

    Here the structural matrix is split whereas the fluid matrix is taken as
    it is.

    \author Shadan Shahmiri
    \date  11/2011
  */
  class FluidFluidMonolithicStructureSplit : public Monolithic
  {
  public:
    explicit FluidFluidMonolithicStructureSplit(Epetra_Comm& comm);

    /*! do the setup for the monolithic system

    <pre>

    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix

    </pre>

    */
    virtual void SetupSystem();

    void SetupNewSystem();

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall);

    /// setup composed system matrix from field solvers
    void SetupSystemMatrixold(LINALG::BlockSparseMatrixBase& mat);

    virtual void SetupSystemMatrix();

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

    virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

  protected:

    /// extract the three field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                     Teuchos::RCP<const Epetra_Vector>& ax);

    virtual void Update();

    virtual void PrepareTimeStep();


  private:

    /// build block vector from field vectors
    void SetupVector(Epetra_Vector &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     Teuchos::RCP<const Epetra_Vector> av,
                     double fluidscale);

    /// block system matrix
    //Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// coupling of fluid and ale (interface only)
    ADAPTER::Coupling icoupfa_;

    /// coupling of fluid and ale at the free surface
    ADAPTER::Coupling fscoupfa_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    UTILS::MatrixRowColTransform sggtransform_;
    UTILS::MatrixRowTransform sgitransform_;
    UTILS::MatrixColTransform sigtransform_;
    UTILS::MatrixColTransform aigtransform_;

    UTILS::MatrixColTransform fmiitransform_;
    UTILS::MatrixColTransform fmgitransform_;

    UTILS::MatrixColTransform fsaigtransform_;
    UTILS::MatrixColTransform fsmgitransform_;

    ///@}

    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> srowsum_;
    Teuchos::RCP<Epetra_Vector> scolsum_;
    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    /// additional ale residual to avoid incremental ale errors
    Teuchos::RCP<Epetra_Vector> aleresidual_;

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;

    FLD::UTILS::FluidXFluidMapExtractor xfluidfluidsplitter_;
  };
}

#endif
#endif
