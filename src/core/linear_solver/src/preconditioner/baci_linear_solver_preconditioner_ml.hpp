/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 0

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_ML_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_ML_HPP

#include "baci_config.hpp"

#include "baci_linear_solver_preconditioner_type.hpp"

BACI_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  /// ml preconditioners
  /*!
    Set of single-matrix algebraic multi-grid preconditioners.
   */
  class MLPreconditioner : public PreconditionerType
  {
   public:
    MLPreconditioner(Teuchos::ParameterList& mllist);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override { return P_; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string getParameterListName() const override { return "ML Parameters"; }

   private:
    Teuchos::ParameterList& mllist_;

    //! system of equations used for preconditioning used by P_ only
    Teuchos::RCP<Epetra_RowMatrix> Pmatrix_;

    /// preconditioner
    Teuchos::RCP<Epetra_Operator> P_;
  };
}  // namespace CORE::LINEAR_SOLVER

BACI_NAMESPACE_CLOSE

#endif
