/*----------------------------------------------------------------------*/
/*! \file

\level 1


\brief Preconditioner for FSI problems with additional constraints
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FSI_CONSTR_OVERLAPPREC_HPP
#define BACI_FSI_CONSTR_OVERLAPPREC_HPP

#include "baci_config.hpp"

#include "baci_fsi_overlapprec.hpp"

BACI_NAMESPACE_OPEN

// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
// #define BLOCKMATRIXMERGE

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
    ConstrOverlappingBlockMatrix(const CORE::LINALG::MultiMapExtractor& maps,
        ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
        bool structuresplit, int symmetric, double omega = 1.0, int iterations = 1,
        double somega = 1.0, int siterations = 0, double fomega = 1.0, int fiterations = 0,
        double aomega = 1.0, int aiterations = 0);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    const char* Label() const override;

    //@}

    /// setup of block preconditioners
    void SetupPreconditioner() override;


   protected:
    /// symmetric Gauss-Seidel block preconditioner
    void SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    //     Teuchos::RCP<CORE::LINALG::SparseMatrix> interconA_;
    //     Teuchos::RCP<Epetra_Vector> interconsol_;
    //     Teuchos::RCP<Epetra_Vector> interconrhs_;
    //     Teuchos::RCP<Epetra_LinearProblem> linprob_
    //     Teuchos::RCP<Amesos_Umfpack> constraintsolver_;

    Teuchos::RCP<Epetra_Map> overallfsimap_;
    CORE::LINALG::MultiMapExtractor fsiextractor_;
    //     Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> invDiag_;

    double alpha_;                 /// "relaxation" parameter in SIMPLE approximation of matrix
    int simpleiter_;               /// number of iterations in SIMPLE preconditioner
    INPAR::FSI::PrecConstr prec_;  /// preconditioner for constraint system
  };
}  // namespace FSI

BACI_NAMESPACE_CLOSE

#endif
