/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_BLOCK_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_BLOCK_HPP

#include "4C_config.hpp"

#include "4C_linear_solver_preconditioner_type.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver
{
  /// SIMPLE(R) block preconditioner
  /*!
    Block preconditioners assume the Epetra_Operator to be a
    Core::LinAlg::BlockSparseMatrix.
   */
  class SimplePreconditioner : public PreconditionerTypeBase
  {
   public:
    SimplePreconditioner(Teuchos::ParameterList& params);

    void setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> prec_operator() const override { return p_; }

   private:
    Teuchos::ParameterList& params_;
    Teuchos::RCP<Epetra_Operator> p_;
  };

  /// General purpose block gauss-seidel preconditioner
  /*!
    2x2 block preconditioner
   */
  class BGSPreconditioner : public PreconditionerTypeBase
  {
   public:
    BGSPreconditioner(Teuchos::ParameterList& params, Teuchos::ParameterList& bgslist);

    void setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> prec_operator() const override { return p_; }

   private:
    Teuchos::ParameterList& params_;
    Teuchos::ParameterList& bgslist_;
    Teuchos::RCP<Epetra_Operator> p_;
  };
}  // namespace Core::LinearSolver

FOUR_C_NAMESPACE_CLOSE

#endif
