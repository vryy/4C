#ifdef CCADISCRET

#ifndef FSI_OVERLAPPREC_H
#define FSI_OVERLAPPREC_H

#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_lib/linalg_precond.H"
#include "../drt_lib/linalg_solver.H"

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

    /// lower Gauss-Seidel block preconditioner
    void LowerGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// upper Gauss-Seidel block preconditioner
    void UpperGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// symmetric Gauss-Seidel block preconditioner
    void SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

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
