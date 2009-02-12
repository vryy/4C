#ifdef CCADISCRET

#ifndef FSI_LAGRANGEPREC_H
#define FSI_LAGRANGEPREC_H

#include "../drt_lib/linalg_blocksparsematrix.H"
#include "../drt_lib/linalg_precond.H"
#include "../drt_lib/linalg_solver.H"

// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
//#define BLOCKMATRIXMERGE

namespace FSI
{

  /// special version of block matrix that includes the FSI block preconditioner
  class LagrangianBlockMatrix : public LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>
  {
  public:

    /// construction
    LagrangianBlockMatrix(const LINALG::MultiMapExtractor& maps,
                          Teuchos::RCP<LINALG::Solver> structuresolver,
                          Teuchos::RCP<LINALG::Solver> fluidsolver,
                          Teuchos::RCP<LINALG::Solver> alesolver);

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

    void SAFLowerGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    // evil debug helper
    void MergeSolve(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /** \name Field solver objects */
    //@{

    Teuchos::RCP<LINALG::Preconditioner> structuresolver_;
    Teuchos::RCP<LINALG::Preconditioner> fluidsolver_;
    Teuchos::RCP<LINALG::Preconditioner> alesolver_;

    //@}
  };

}

#endif
#endif
