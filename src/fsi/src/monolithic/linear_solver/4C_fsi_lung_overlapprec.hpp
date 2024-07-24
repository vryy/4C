/*----------------------------------------------------------------------*/
/*! \file
\brief BGS preconditioner for volume-coupled FSI
\level 2
*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FSI_LUNG_OVERLAPPREC_HPP
#define FOUR_C_FSI_LUNG_OVERLAPPREC_HPP

#include "4C_config.hpp"

#include "4C_fsi_overlapprec.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
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
    Teuchos::RCP<Core::LinAlg::SparseMatrix> calculate_schur(const Core::LinAlg::SparseMatrix& A,
        const Core::LinAlg::SparseMatrix& B, const Core::LinAlg::SparseMatrix& C);

   private:
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> res_;
  };


  /// special version of block matrix that includes the FSI block
  /// preconditioner as well as a SIMPLE preconditioner for handling
  /// the constraint part for lung fsi simulations
  class LungOverlappingBlockMatrix : public OverlappingBlockMatrix
  {
   public:
    /// construction
    LungOverlappingBlockMatrix(const Core::LinAlg::MultiMapExtractor& maps,
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

    Teuchos::RCP<LungSchurComplement> StructSchur_;
    Teuchos::RCP<LungSchurComplement> FluidSchur_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> interconA_;
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> invDiag_;

    Teuchos::RCP<Core::LinAlg::Solver> constraintsolver_;
    Teuchos::RCP<Epetra_Map> overallfsimap_;
    Core::LinAlg::MultiMapExtractor fsiextractor_;

    double alpha_;                 /// "relaxation" parameter in SIMPLE approximation of matrix
    int simpleiter_;               /// number of iterations in SIMPLE preconditioner
    Inpar::FSI::PrecConstr prec_;  /// preconditioner for constraint system
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
