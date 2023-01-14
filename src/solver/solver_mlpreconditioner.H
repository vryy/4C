/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 0

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_MLPRECONDITIONER_H_
#define SOLVER_MLPRECONDITIONER_H_

#include "solver_preconditionertype.H"

namespace LINALG::SOLVER
{
  /// ml preconditioners
  /*!
    Set of single-matrix algebraic multi-grid preconditioners.
   */
  class MLPreconditioner : public PreconditionerType
  {
   public:
    MLPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    Epetra_Operator* PrecOperator() const override { return &*P_; }

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override { return P_; }

    void Print(std::ostream& stream) const override { stream << "ML"; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    const std::string getParameterListName() const override { return "ML Parameters"; }

   private:
    Teuchos::ParameterList& mllist_;

    //! system of equations used for preconditioning used by P_ only
    Teuchos::RCP<Epetra_RowMatrix> Pmatrix_;

    /// preconditioner
    Teuchos::RCP<Epetra_Operator> P_;
  };
}  // namespace LINALG::SOLVER

#endif
