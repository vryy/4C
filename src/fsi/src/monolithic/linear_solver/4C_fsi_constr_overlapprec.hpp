/*----------------------------------------------------------------------*/
/*! \file

\level 1


\brief Preconditioner for FSI problems with additional constraints
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_CONSTR_OVERLAPPREC_HPP
#define FOUR_C_FSI_CONSTR_OVERLAPPREC_HPP

#include "4C_config.hpp"

#include "4C_fsi_overlapprec.hpp"

FOUR_C_NAMESPACE_OPEN

// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
// #define BLOCKMATRIXMERGE

namespace Adapter
{
  class Structure;
  class Fluid;
}  // namespace Adapter

namespace FSI
{
  /// special version of block matrix that includes the FSI block
  /// preconditioner as well as a SIMPLE preconditioner for handling
  /// the constraint part for lung fsi simulations
  class ConstrOverlappingBlockMatrix : public OverlappingBlockMatrix
  {
   public:
    /// construction
    ConstrOverlappingBlockMatrix(const Core::LinAlg::MultiMapExtractor& maps,
        Adapter::FSIStructureWrapper& structure, Adapter::Fluid& fluid, Adapter::AleFsiWrapper& ale,
        bool structuresplit, int symmetric, double omega = 1.0, int iterations = 1,
        double somega = 1.0, int siterations = 0, double fomega = 1.0, int fiterations = 0,
        double aomega = 1.0, int aiterations = 0);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    const char* Label() const override;

    //@}

    /// setup of block preconditioners
    void setup_preconditioner() override;


   protected:
    /// symmetric Gauss-Seidel block preconditioner
    void sgs(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    //     Teuchos::RCP<Core::LinAlg::SparseMatrix> interconA_;
    //     Teuchos::RCP<Epetra_Vector> interconsol_;
    //     Teuchos::RCP<Epetra_Vector> interconrhs_;
    //     Teuchos::RCP<Epetra_LinearProblem> linprob_
    //     Teuchos::RCP<Amesos_Umfpack> constraintsolver_;

    Teuchos::RCP<Epetra_Map> overallfsimap_;
    Core::LinAlg::MultiMapExtractor fsiextractor_;
    //     Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> invDiag_;

    double alpha_;                 /// "relaxation" parameter in SIMPLE approximation of matrix
    int simpleiter_;               /// number of iterations in SIMPLE preconditioner
    Inpar::FSI::PrecConstr prec_;  /// preconditioner for constraint system
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
