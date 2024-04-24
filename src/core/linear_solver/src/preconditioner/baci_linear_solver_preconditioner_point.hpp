/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_POINT_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_POINT_HPP

#include "baci_config.hpp"

#include "baci_linalg_downwindmatrix.hpp"
#include "baci_linear_solver_preconditioner_type.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  /// inf-norm scaling fits the preconditioner framework perfecty
  /*!
   Modifies the underlying matrix and needs to unscale the result
   afterwards. Can be combined with any single-matrix preconditioner.
  */
  class InfNormPreconditioner : public PreconditionerTypeBase
  {
   public:
    InfNormPreconditioner(Teuchos::RCP<PreconditionerTypeBase> preconditioner);

    Epetra_LinearProblem& LinearProblem() override { return preconditioner_->LinearProblem(); }

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override
    {
      return preconditioner_->PrecOperator();
    }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string getParameterListName() const override { return "unknown"; }

   private:
    Teuchos::RCP<PreconditionerTypeBase> preconditioner_;
    Teuchos::RCP<Epetra_Vector> rowsum_;
    Teuchos::RCP<Epetra_Vector> colsum_;
  };

  /// diagonal scaling fits the preconditioner framework perfecty
  /*!
   Modifies the underlying matrix and needs to unscale the result
   afterwards. Can be combined with any single-matrix preconditioner.
  */
  class SymDiagPreconditioner : public PreconditionerTypeBase
  {
   public:
    SymDiagPreconditioner(Teuchos::RCP<PreconditionerTypeBase> preconditioner);

    Epetra_LinearProblem& LinearProblem() override { return preconditioner_->LinearProblem(); }

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override
    {
      return preconditioner_->PrecOperator();
    }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string getParameterListName() const override { return "unknown"; }

   private:
    /// embedded preconditioner
    Teuchos::RCP<LINEAR_SOLVER::PreconditionerTypeBase> preconditioner_;

    /// matrix diagonal
    Teuchos::RCP<Epetra_Vector> diag_;
  };
}  // namespace CORE::LINEAR_SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
