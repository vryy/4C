#ifdef CCADISCRET

#ifndef FSI_OVERLAPPREC_H
#define FSI_OVERLAPPREC_H

#include "../drt_lib/linalg_blocksparsematrix.H"
#include "../drt_lib/linalg_precond.H"
#include "../drt_lib/linalg_solver.H"

#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_fluid.H"
#include "../drt_adapter/adapter_ale.H"

// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
#define BLOCKMATRIXMERGE

namespace FSI
{

  /// Base class for all FSI block preconditioning matrices
  class BlockPreconditioningMatrix : public LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>
  {
  public:

    BlockPreconditioningMatrix(const LINALG::MultiMapExtractor& maps,
                               ADAPTER::Structure& structure,
                               ADAPTER::Fluid& fluid,
                               ADAPTER::Ale& ale,
                               int symmetric,
                               double omega=1.0,
                               int iterations=1,
                               double somega=1.0,
                               int siterations=0,
                               double fomega=1.0,
                               int fiterations=0,
                               FILE* err=NULL);

    /** \name Mathematical functions */
    //@{

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    virtual int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    //@}

    /// setup of block preconditioners
    virtual void SetupPreconditioner();

  protected:

    /// (symmetric) Gauss-Seidel block preconditioner
    virtual void SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const = 0;

    /// merge block matrix for direct solve
    void MergeSolve(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// Richardson iteration on one block using the given flags
    static void LocalBlockRichardson(Teuchos::RCP<LINALG::Preconditioner> solver,
                                     const LINALG::SparseMatrix& innerOp,
                                     Teuchos::RCP<Epetra_Vector> x,
                                     Teuchos::RCP<Epetra_Vector> y,
                                     Teuchos::RCP<Epetra_Vector> tmpx,
                                     int iterations,
                                     double omega,
                                     FILE* err,
                                     const Epetra_Comm& comm);

    /** \name Field solver objects */
    //@{

    Teuchos::RCP<LINALG::Preconditioner> structuresolver_;
    Teuchos::RCP<LINALG::Preconditioner> fluidsolver_;
    Teuchos::RCP<LINALG::Preconditioner> alesolver_;

    Teuchos::RCP<LINALG::Preconditioner> constalesolver_;

    //@}

    /// Symmetric block GS preconditioner in monolithic FSI or ordinary GS
    int symmetric_;

    /// \name Richardson iteration
    //@{

    double omega_;
    int iterations_;
    double somega_;
    int siterations_;
    double fomega_;
    int fiterations_;

    //@}

    /// log file
    FILE* err_;

#ifdef BLOCKMATRIXMERGE
    /// debug merged sparse
    Teuchos::RCP<LINALG::SparseMatrix> sparse_;
#endif
  };


  /// special version of block matrix that includes the FSI block preconditioner
  /*!
    The normal block matrix is enhanced by a ApplyInverse() method that does
    the Gauss-Seidel block preconditioning explicitly for FSI block matrices.
   */
  class OverlappingBlockMatrix : public BlockPreconditioningMatrix
  {
  public:

    /// construction
    OverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
                           ADAPTER::Structure& structure,
                           ADAPTER::Fluid& fluid,
                           ADAPTER::Ale& ale,
                           bool structuresplit,
                           int symmetric,
                           double omega=1.0,
                           int iterations=1,
                           double somega=1.0,
                           int siterations=0,
                           double fomega=1.0,
                           int fiterations=0,
                           FILE* err=NULL);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

  private:

    /// symmetric Gauss-Seidel block preconditioner
    void SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// split is in structural matrix, interface equations belong to fluid block
    bool structuresplit_;
  };
}

#endif
#endif
