/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_POINTPRECONDITIONER_H_
#define SOLVER_POINTPRECONDITIONER_H_

#include "../linalg/linalg_downwindmatrix.H"
#include "../linalg/linalg_none_operator.H"

#include "solver_preconditionertype.H"

namespace LINALG::SOLVER
{
  class NonePreconditioner : public PreconditionerType
  {
   public:
    NonePreconditioner(FILE* outfile, Teuchos::ParameterList& list);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    Epetra_Operator* PrecOperator() const override { return &*prec_; }

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override { return prec_; }

    void Print(std::ostream& stream) override { stream << "None"; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    const std::string getParameterListName() const override { return "unknown"; }

   private:
    Teuchos::RCP<LINALG::NoneOperator> prec_;
  };

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

    void Print(std::ostream& stream) override
    {
      stream << "InfNorm(";
      preconditioner_->Print(stream);
      stream << ")";
    }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    const std::string getParameterListName() const override { return "unknown"; }

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

    void Print(std::ostream& stream) override
    {
      stream << "SymDiag(";
      preconditioner_->Print(stream);
      stream << ")";
    }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    const std::string getParameterListName() const override { return "unknown"; }

   private:
    /// embedded preconditioner
    Teuchos::RCP<SOLVER::PreconditionerType> preconditioner_;

    /// matrix diagonal
    Teuchos::RCP<Epetra_Vector> diag_;
  };
}  // namespace LINALG::SOLVER

#endif
