/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_BGS2X2_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_BGS2X2_HPP

#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linear_solver_preconditioner_linalg.hpp"
#include "4C_linear_solver_preconditioner_type.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /// Block Gauss-Seidel preconditioner for a 2x2 system
  class BgS2x2Operator : public Epetra_Operator
  {
   public:
    /*!
    \brief Standard Constructor
    */
    explicit BgS2x2Operator(Teuchos::RCP<Epetra_Operator> A, const Teuchos::ParameterList& list1,
        const Teuchos::ParameterList& list2, int global_iter, double global_omega, int block1_iter,
        double block1_omega, int block2_iter, double block2_omega, bool fliporder);

    /*!
    \brief Destructor
    */
    ~BgS2x2Operator() override
    {
      solver1_ = Teuchos::null;
      solver2_ = Teuchos::null;
    }

    /// Llabel of this class.
    const char* Label() const override { return "Core::LinAlg::BGS2x2_Operator"; }

    /// Comm of this class
    const Epetra_Comm& Comm() const override { return (a_->Comm()); }


    /// Operator domain map
    const Epetra_Map& OperatorDomainMap() const override { return a_->FullDomainMap(); }

    /// Operator range map
    const Epetra_Map& OperatorRangeMap() const override { return a_->FullRangeMap(); }

    /// Setup of preconditioners for individual blocks
    void setup_block_preconditioners();

    /// Apply inverse of the preconditioner
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /// not implemented
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
    {
      FOUR_C_THROW("Apply does not make sense for Core::LinAlg::BGS2x2_Operator");
      return (-1);
    }

    int SetUseTranspose(bool UseTranspose) override
    {
      // we default to false
      return 0;
    }

    /// not implemented
    double NormInf() const override
    {
      FOUR_C_THROW("NormInf not impl.");
      return (-1.0);
    }

    bool UseTranspose() const override
    {
      // we default to false
      return false;
    }

    /// not implemented
    bool HasNormInf() const override
    {
      FOUR_C_THROW("HasNormInf not impl.");
      return false;
    }

   private:
    // don't want copy-ctor and = operator
    BgS2x2Operator(BgS2x2Operator& old);
    BgS2x2Operator operator=(const BgS2x2Operator& old);

    /// Richardson iteration on one block using the given flags
    void local_block_richardson(Teuchos::RCP<Preconditioner> solver, const SparseMatrix& Op,
        Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> y,
        Teuchos::RCP<Epetra_MultiVector> tmpx, int iter, double omega) const;


    Teuchos::ParameterList list1_;  // list for solver of first diagonal block
    Teuchos::ParameterList list2_;  // list for solver of second diagonal block

    MultiMapExtractor mmex_;                 // a  multimapetxractor to handle extracts
    Teuchos::RCP<BlockSparseMatrixBase> a_;  // 2x2 block matrix

    Teuchos::RCP<Preconditioner> solver1_;  // solver of block 1
    Teuchos::RCP<Preconditioner> solver2_;  // solver of block 2

    int global_iter_;
    double global_omega_;
    int block1_iter_;
    double block1_omega_;
    int block2_iter_;
    double block2_omega_;

    int firstind_;  /// index of block "0" in Gauss-Seidel procedure
    int secind_;    /// index of block "1" in Gauss-Seidel procedure
  };
}  // namespace Core::LinAlg


FOUR_C_NAMESPACE_CLOSE

#endif
