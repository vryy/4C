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
                           Teuchos::RCP<LINALG::Solver> alesolver,
                           bool structuresplit=false,
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

    /// fluid-structure-ale lower Gauss-Seidel block preconditioner
    void FSALowerGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// upper Gauss-Seidel block preconditioner
    void UpperGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// symmetric Gauss-Seidel block preconditioner
    void SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// merge block matrix for direct solve
    void MergeSolve(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// block Richardson iteration with fixed relaxation
    void BlockRichardson(const LINALG::SparseMatrix& Aii,
                         const LINALG::SparseMatrix& Aig,
                         const LINALG::SparseMatrix& Agi,
                         const LINALG::SparseMatrix& Agg,
                         Teuchos::RCP<LINALG::Preconditioner> Sii,
                         Teuchos::RCP<LINALG::Preconditioner> Sgg,
                         Teuchos::RCP<Epetra_Vector> xi,
                         Teuchos::RCP<Epetra_Vector> xg,
                         Teuchos::RCP<Epetra_Vector> bi,
                         Teuchos::RCP<Epetra_Vector> bg,
                         double omega,
                         int itenum) const;

    /** \name Field solver objects */
    //@{

    Teuchos::RCP<LINALG::Preconditioner> structuresolver_;
    Teuchos::RCP<LINALG::Preconditioner> fluidsolver_;
    Teuchos::RCP<LINALG::Preconditioner> alesolver_;

    //@}

    /// split is in structural matrix, interface equations belong to fluid block
    bool structuresplit_;

    /// \name Richardson iteration
    //@{

    double somega_;
    int siterations_;
    double fomega_;
    int fiterations_;

    //@}

    /// log file
    FILE* err_;
  };
}

#endif
#endif
