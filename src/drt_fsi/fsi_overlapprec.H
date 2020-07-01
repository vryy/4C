/*----------------------------------------------------------------------*/
/*! \file

\brief Base class for all FSI block preconditioning matrices

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FSI_OVERLAPPREC_H
#define FSI_OVERLAPPREC_H

#include "../linalg/linalg_blocksparsematrix.H"

#include "../drt_inpar/inpar_fsi.H"


// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
//#define BLOCKMATRIXMERGE

// forward declarations
namespace ADAPTER
{
  class AleFsiWrapper;
  class Fluid;
  class FSIStructureWrapper;
}  // namespace ADAPTER

namespace FSI
{
  namespace UTILS
  {
    class MonolithicDebugWriter;
  }
}  // namespace FSI

namespace LINALG
{
  class Preconditioner;
}

namespace FSI
{
  /// Base class for all FSI block preconditioning matrices
  class BlockPreconditioningMatrix
      : public LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>
  {
   public:
    BlockPreconditioningMatrix(Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg,
        const LINALG::MultiMapExtractor& maps, ::ADAPTER::FSIStructureWrapper& structure,
        ::ADAPTER::Fluid& fluid, ::ADAPTER::AleFsiWrapper& ale, int symmetric, double omega = 1.0,
        int iterations = 1, double somega = 1.0, int siterations = 0, double fomega = 1.0,
        int fiterations = 0, double aomega = 1.0, int aiterations = 0, FILE* err = NULL);

    /** \name Mathematical functions */
    //@{

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //@}

    /// setup of block preconditioners
    virtual void SetupPreconditioner();

   protected:
    /// (symmetric) Gauss-Seidel block preconditioner
    virtual void SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

    /// merge block matrix for direct solve
    void MergeSolve(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// Richardson iteration on one block using the given flags
    static void LocalBlockRichardson(Teuchos::RCP<LINALG::Preconditioner> solver,
        const LINALG::SparseMatrix& innerOp, Teuchos::RCP<Epetra_Vector> x,
        Teuchos::RCP<Epetra_Vector> y, Teuchos::RCP<Epetra_Vector> tmpx, int iterations,
        double omega, FILE* err, const Epetra_Comm& comm);

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
    double aomega_;
    int aiterations_;

    //@}

    /// log file
    FILE* err_;

    /// debug writer
    Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg_;

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
    OverlappingBlockMatrix(Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg,
        const LINALG::MultiMapExtractor& maps, ::ADAPTER::FSIStructureWrapper& structure,
        ::ADAPTER::Fluid& fluid, ::ADAPTER::AleFsiWrapper& ale, bool structuresplit, int symmetric,
        double omega = 1.0, int iterations = 1, double somega = 1.0, int siterations = 0,
        double fomega = 1.0, int fiterations = 0, double aomega = 1.0, int aiterations = 0,
        FILE* err = NULL);

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

    /// split is in structural matrix, interface equations belong to fluid block
    bool structuresplit_;

    ::ADAPTER::FSIStructureWrapper& structure_;
    ::ADAPTER::Fluid& fluid_;
    ::ADAPTER::AleFsiWrapper& ale_;
  };
}  // namespace FSI

#endif
