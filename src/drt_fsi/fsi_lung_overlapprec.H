/*----------------------------------------------------------------------*/
/*!
\brief BGS preconditioner for volume-coupled FSI
\level 2
\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_LUNG_OVERLAPPREC_H
#define FSI_LUNG_OVERLAPPREC_H

#include "fsi_overlapprec.H"

// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
//#define BLOCKMATRIXMERGE

// forward declarations
namespace LINALG
{
  class Solver;
}

namespace FSI
{
  /// this helper class is needed to save the graph of a temporary matrix and the
  /// Schur complement -> the method "CalculateSchur" needs to be called
  /// always with the same three matrices!
  class LungSchurComplement
  {
   public:
    /// construction
    LungSchurComplement(){};

    /// determination of the Schur complement
    Teuchos::RCP<LINALG::SparseMatrix> CalculateSchur(const LINALG::SparseMatrix& A,
        const LINALG::SparseMatrix& B, const LINALG::SparseMatrix& C);

   private:
    Teuchos::RCP<LINALG::SparseMatrix> temp_;
    Teuchos::RCP<LINALG::SparseMatrix> res_;
  };


  /// special version of block matrix that includes the FSI block
  /// preconditioner as well as a SIMPLE preconditioner for handling
  /// the constraint part for lung fsi simulations
  class LungOverlappingBlockMatrix : public OverlappingBlockMatrix
  {
   public:
    /// construction
    LungOverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
        ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
        bool structuresplit, int symmetric, double omega = 1.0, int iterations = 1,
        double somega = 1.0, int siterations = 0, double fomega = 1.0, int fiterations = 0,
        double aomega = 1.0, int aiterations = 0, FILE* err = NULL);

    /// destruction
    ~LungOverlappingBlockMatrix(){};

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

    Teuchos::RCP<LungSchurComplement> StructSchur_;
    Teuchos::RCP<LungSchurComplement> FluidSchur_;
    Teuchos::RCP<LINALG::SparseMatrix> interconA_;
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> invDiag_;

    Teuchos::RCP<LINALG::Solver> constraintsolver_;
    Teuchos::RCP<Epetra_Map> overallfsimap_;
    LINALG::MultiMapExtractor fsiextractor_;

    double alpha_;                 /// "relaxation" parameter in SIMPLE approximation of matrix
    int simpleiter_;               /// number of iterations in SIMPLE preconditioner
    INPAR::FSI::PrecConstr prec_;  /// preconditioner for constraint system
  };
}  // namespace FSI

#endif
