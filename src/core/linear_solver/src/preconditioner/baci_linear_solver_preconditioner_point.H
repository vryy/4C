/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_PRECONDITIONER_POINT_H
#define BACI_LINEAR_SOLVER_PRECONDITIONER_POINT_H

#include "baci_linalg_downwindmatrix.H"
#include "baci_linear_solver_preconditioner_type.H"

namespace CORE::LINEAR_SOLVER
{
  /// inf-norm scaling fits the preconditioner framework perfecty
  /*!
   Modifies the underlying matrix and needs to unscale the result
   afterwards. Can be combined with any single-matrix preconditioner.
  */
  class InfNormPreconditioner : public PreconditionerType
  {
   public:
    InfNormPreconditioner(Teuchos::RCP<PreconditionerType> preconditioner);

    Epetra_LinearProblem& LinearProblem() override { return preconditioner_->LinearProblem(); }

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) override;

    Epetra_Operator* PrecOperator() const override { return preconditioner_->PrecOperator(); }

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override
    {
      return preconditioner_->PrecOperatorRCP();
    }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string getParameterListName() const override { return "unknown"; }

   private:
    Teuchos::RCP<PreconditionerType> preconditioner_;
    Teuchos::RCP<Epetra_Vector> rowsum_;
    Teuchos::RCP<Epetra_Vector> colsum_;
  };

  /// diagonal scaling fits the preconditioner framework perfecty
  /*!
   Modifies the underlying matrix and needs to unscale the result
   afterwards. Can be combined with any single-matrix preconditioner.
  */
  class SymDiagPreconditioner : public PreconditionerType
  {
   public:
    SymDiagPreconditioner(Teuchos::RCP<PreconditionerType> preconditioner);

    Epetra_LinearProblem& LinearProblem() override { return preconditioner_->LinearProblem(); }

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) override;

    Epetra_Operator* PrecOperator() const override { return preconditioner_->PrecOperator(); }

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override
    {
      return preconditioner_->PrecOperatorRCP();
    }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string getParameterListName() const override { return "unknown"; }

   private:
    /// embedded preconditioner
    Teuchos::RCP<LINEAR_SOLVER::PreconditionerType> preconditioner_;

    /// matrix diagonal
    Teuchos::RCP<Epetra_Vector> diag_;
  };
}  // namespace CORE::LINEAR_SOLVER

#endif
