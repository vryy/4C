/*----------------------------------------------------------------------*/
/*! \file

\brief Base class for all FSI block preconditioning matrices

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_OVERLAPPREC_HPP
#define FOUR_C_FSI_OVERLAPPREC_HPP

#include "4C_config.hpp"

#include "4C_inpar_fsi.hpp"
#include "4C_linalg_blocksparsematrix.hpp"


FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Adapter
{
  class AleFsiWrapper;
  class Fluid;
  class FSIStructureWrapper;
}  // namespace Adapter

namespace FSI
{
  namespace UTILS
  {
    class MonolithicDebugWriter;
  }
}  // namespace FSI

namespace Core::LinAlg
{
  class Preconditioner;
}

namespace FSI
{
  /// Base class for all FSI block preconditioning matrices
  class BlockPreconditioningMatrix
      : public Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>
  {
   public:
    BlockPreconditioningMatrix(Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg,
        const Core::LinAlg::MultiMapExtractor& maps, Adapter::FSIStructureWrapper& structure,
        Adapter::Fluid& fluid, Adapter::AleFsiWrapper& ale, int symmetric, double omega = 1.0,
        int iterations = 1, double somega = 1.0, int siterations = 0, double fomega = 1.0,
        int fiterations = 0, double aomega = 1.0, int aiterations = 0, FILE* err = nullptr);

    /** \name Mathematical functions */
    //@{

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    //@}

   protected:
    /// (symmetric) Gauss-Seidel block preconditioner
    virtual void sgs(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

    /// Richardson iteration on one block using the given flags
    static void local_block_richardson(Teuchos::RCP<Core::LinAlg::Preconditioner> solver,
        const Core::LinAlg::SparseMatrix& innerOp, Teuchos::RCP<Epetra_Vector> x,
        Teuchos::RCP<Epetra_Vector> y, Teuchos::RCP<Epetra_Vector> tmpx, int iterations,
        double omega, FILE* err, const Epetra_Comm& comm);

    /** \name Field solver objects */
    //@{

    Teuchos::RCP<Core::LinAlg::Preconditioner> structuresolver_;
    Teuchos::RCP<Core::LinAlg::Preconditioner> fluidsolver_;
    Teuchos::RCP<Core::LinAlg::Preconditioner> alesolver_;

    Teuchos::RCP<Core::LinAlg::Preconditioner> constalesolver_;

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
        const Core::LinAlg::MultiMapExtractor& maps, Adapter::FSIStructureWrapper& structure,
        Adapter::Fluid& fluid, Adapter::AleFsiWrapper& ale, bool structuresplit, int symmetric,
        double omega = 1.0, int iterations = 1, double somega = 1.0, int siterations = 0,
        double fomega = 1.0, int fiterations = 0, double aomega = 1.0, int aiterations = 0);

   protected:
    /// symmetric Gauss-Seidel block preconditioner
    void sgs(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override{};

    /// split is in structural matrix, interface equations belong to fluid block
    bool structuresplit_;

    Adapter::FSIStructureWrapper& structure_;
    Adapter::Fluid& fluid_;
    Adapter::AleFsiWrapper& ale_;
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
