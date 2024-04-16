/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_BGS2X2_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_BGS2X2_HPP

#include "baci_config.hpp"

#include "baci_linalg_blocksparsematrix.hpp"
#include "baci_linalg_mapextractor.hpp"
#include "baci_linear_solver_preconditioner_linalg.hpp"
#include "baci_linear_solver_preconditioner_type.hpp"
#include "baci_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  /// Block Gauss-Seidel preconditioner for a 2x2 system
  class BGS2x2_Operator : public Epetra_Operator
  {
   public:
    /*!
    \brief Standard Constructor
    */
    explicit BGS2x2_Operator(Teuchos::RCP<Epetra_Operator> A, const Teuchos::ParameterList& list1,
        const Teuchos::ParameterList& list2, int global_iter, double global_omega, int block1_iter,
        double block1_omega, int block2_iter, double block2_omega, bool fliporder);

    /*!
    \brief Destructor
    */
    ~BGS2x2_Operator() override
    {
      solver1_ = Teuchos::null;
      solver2_ = Teuchos::null;
    }

    /// Llabel of this class.
    const char* Label() const override { return "CORE::LINALG::BGS2x2_Operator"; }

    /// Comm of this class
    const Epetra_Comm& Comm() const override { return (A_->Comm()); }


    /// Operator domain map
    const Epetra_Map& OperatorDomainMap() const override { return A_->FullDomainMap(); }

    /// Operator range map
    const Epetra_Map& OperatorRangeMap() const override { return A_->FullRangeMap(); }

    /// Setup of preconditioners for individual blocks
    void SetupBlockPreconditioners();

    /// Apply inverse of the preconditioner
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /// not implemented
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
    {
      dserror("Apply does not make sense for CORE::LINALG::BGS2x2_Operator");
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
      dserror("NormInf not impl.");
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
      dserror("HasNormInf not impl.");
      return false;
    }

   private:
    // don't want copy-ctor and = operator
    BGS2x2_Operator(BGS2x2_Operator& old);
    BGS2x2_Operator operator=(const BGS2x2_Operator& old);

    /// Richardson iteration on one block using the given flags
    void LocalBlockRichardson(Teuchos::RCP<Preconditioner> solver, const SparseMatrix& Op,
        Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> y,
        Teuchos::RCP<Epetra_MultiVector> tmpx, int iter, double omega) const;


    Teuchos::ParameterList list1_;  // list for solver of first diagonal block
    Teuchos::ParameterList list2_;  // list for solver of second diagonal block

    MultiMapExtractor mmex_;                 // a  multimapetxractor to handle extracts
    Teuchos::RCP<BlockSparseMatrixBase> A_;  // 2x2 block matrix

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
}  // namespace CORE::LINALG


BACI_NAMESPACE_CLOSE

#endif
