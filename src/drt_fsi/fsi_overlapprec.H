#ifdef CCADISCRET

#ifndef FSI_OVERLAPPREC_H
#define FSI_OVERLAPPREC_H

#include "../drt_lib/linalg_systemmatrix.H"
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

  private:

    /** \name Field solver objects */
    //@{

    Teuchos::RCP<LINALG::Solver> structuresolver_;
    Teuchos::RCP<LINALG::Solver> fluidsolver_;
    Teuchos::RCP<LINALG::Solver> alesolver_;

    //@}
  };
}

#endif
#endif
