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
//#define BLOCKMATRIXMERGE

namespace FSI
{

  /// special version of block matrix that includes the FSI block preconditioner
  /*!
    The normal block matrix is enhanced by a ApplyInverse() method that does
    the Gauss-Seidel block preconditioning explicitly for FSI block matrices.
   */
  class OverlappingBlockMatrix : public LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>
  {
  public:

    /// construction
    OverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
                           ADAPTER::Structure& structure,
                           ADAPTER::Fluid& fluid,
                           ADAPTER::Ale& ale,
                           bool structuresplit,
                           int symmetric,
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

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

    /// setup of block preconditioners
    void SetupPreconditioner();

  private:

    /// structure-ale-fluid lower Gauss-Seidel block preconditioner
    void SAFLowerGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// symmetric Gauss-Seidel block preconditioner
    void SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// merge block matrix for direct solve
    void MergeSolve(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /** \name Field solver objects */
    //@{

    Teuchos::RCP<LINALG::Preconditioner> structuresolver_;
    Teuchos::RCP<LINALG::Preconditioner> fluidsolver_;
    Teuchos::RCP<LINALG::Preconditioner> alesolver_;

    Teuchos::RCP<LINALG::Preconditioner> constalesolver_;

    //@}

    /// split is in structural matrix, interface equations belong to fluid block
    bool structuresplit_;

    /// Symmetric block GS preconditioner in monolithic FSI or ordinary GS
    int symmetric_;

    /// \name Richardson iteration
    //@{

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
}

#endif
#endif
