

#ifndef FSI_CONSTR_OVERLAPPREC_H
#define FSI_CONSTR_OVERLAPPREC_H

#include "fsi_overlapprec.H"

// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
//#define BLOCKMATRIXMERGE

namespace ADAPTER
{
  class Structure;
  class Fluid;
}  // namespace ADAPTER

namespace FSI
{
  /// special version of block matrix that includes the FSI block
  /// preconditioner as well as a SIMPLE preconditioner for handling
  /// the constraint part for lung fsi simulations
  class ConstrOverlappingBlockMatrix : public OverlappingBlockMatrix
  {
   public:
    /// construction
    ConstrOverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
        ::ADAPTER::FSIStructureWrapper& structure, ::ADAPTER::Fluid& fluid,
        ::ADAPTER::AleFsiWrapper& ale, bool structuresplit, int symmetric, double omega = 1.0,
        int iterations = 1, double somega = 1.0, int siterations = 0, double fomega = 1.0,
        int fiterations = 0, double aomega = 1.0, int aiterations = 0, FILE* err = NULL);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

    /// setup of block preconditioners
    virtual void SetupPreconditioner();


   protected:
    /// symmetric Gauss-Seidel block preconditioner
    virtual void SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //     Teuchos::RCP<LINALG::SparseMatrix> interconA_;
    //     Teuchos::RCP<Epetra_Vector> interconsol_;
    //     Teuchos::RCP<Epetra_Vector> interconrhs_;
    //     Teuchos::RCP<Epetra_LinearProblem> linprob_
    //     Teuchos::RCP<Amesos_Umfpack> constraintsolver_;

    Teuchos::RCP<Epetra_Map> overallfsimap_;
    LINALG::MultiMapExtractor fsiextractor_;
    //     Teuchos::RCP<LINALG::BlockSparseMatrixBase> invDiag_;

    double alpha_;                 /// "relaxation" parameter in SIMPLE approximation of matrix
    int simpleiter_;               /// number of iterations in SIMPLE preconditioner
    INPAR::FSI::PrecConstr prec_;  /// preconditioner for constraint system
  };
}  // namespace FSI

#endif
