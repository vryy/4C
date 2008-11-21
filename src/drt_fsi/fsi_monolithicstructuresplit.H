#ifdef CCADISCRET

#ifndef FSI_MONOLITHICSTRUCTURESPLIT_H
#define FSI_MONOLITHICSTRUCTURESPLIT_H

#include "../drt_inpar/inpar_fsi.H"
#include "fsi_monolithic.H"
#include "fsi_overlapprec.H"
#include "fsi_matrixtransform.H"

namespace FSI
{

  /// monolithic FSI algorithm with overlapping interface equations
  /*!

    Here the structural matrix is split whereas the fluid matrix is taken as
    it is.

    \sa MonolithicOverlap
    \author u.kue
    \date 06/08
   */
  class MonolithicStructureSplit : public BlockMonolithic
  {
  public:
    explicit MonolithicStructureSplit(Epetra_Comm& comm);

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall=false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x);

  protected:

    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils);

    /// setup of NOX convergance tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Epetra::Group> grp);

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


  private:

    /// build block vector from field vectors
    void SetupVector(Epetra_Vector &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     Teuchos::RCP<const Epetra_Vector> av,
                     double fluidscale);

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// coupling of fluid and ale (interface only)
    ADAPTER::Coupling icoupfa_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    UTILS::MatrixRowColTransform sggtransform_;
    UTILS::MatrixRowTransform sgitransform_;
    UTILS::MatrixColTransform sigtransform_;
    UTILS::MatrixColTransform aigtransform_;

    UTILS::MatrixColTransform fmiitransform_;
    UTILS::MatrixColTransform fmgitransform_;

    ///@}

    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> srowsum_;
    Teuchos::RCP<Epetra_Vector> scolsum_;
    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;
  };
}

#endif
#endif
