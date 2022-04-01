/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 1
Created on: Jul 4, 2011
*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_POINTPRECONDITIONER_H_
#define SOLVER_POINTPRECONDITIONER_H_

#include "../linalg/linalg_downwindmatrix.H"
#include "../linalg/linalg_none_operator.H"

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    class NonePreconditioner : public PreconditionerType
    {
     public:
      NonePreconditioner(FILE* outfile, Teuchos::ParameterList& list);

      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual Epetra_Operator* PrecOperator() const { return &*prec_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return prec_; }

      virtual void Print(std::ostream& stream) { stream << "None"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "unknown"; }

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

      virtual Epetra_LinearProblem& LinearProblem() { return preconditioner_->LinearProblem(); }

      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual Epetra_Operator* PrecOperator() const { return preconditioner_->PrecOperator(); }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const
      {
        return preconditioner_->PrecOperatorRCP();
      }

      virtual void Print(std::ostream& stream)
      {
        stream << "InfNorm(";
        preconditioner_->Print(stream);
        stream << ")";
      }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "unknown"; }

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

      virtual Epetra_LinearProblem& LinearProblem() { return preconditioner_->LinearProblem(); }

      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual Epetra_Operator* PrecOperator() const { return preconditioner_->PrecOperator(); }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const
      {
        return preconditioner_->PrecOperatorRCP();
      }

      virtual void Print(std::ostream& stream)
      {
        stream << "SymDiag(";
        preconditioner_->Print(stream);
        stream << ")";
      }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "unknown"; }

     private:
      /// embedded preconditioner
      Teuchos::RCP<SOLVER::PreconditionerType> preconditioner_;

      /// matrix diagonal
      Teuchos::RCP<Epetra_Vector> diag_;
    };

  }  // namespace SOLVER
}  // namespace LINALG


#endif /* SOLVER_POINTPRECONDITIONER_H_ */
