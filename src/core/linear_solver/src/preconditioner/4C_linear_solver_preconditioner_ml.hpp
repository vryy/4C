/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 0

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_ML_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_ML_HPP

#include "4C_config.hpp"

#include "4C_linear_solver_preconditioner_type.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  /// ml preconditioners
  /*!
    Set of single-matrix algebraic multi-grid preconditioners.
   */
  class MLPreconditioner : public PreconditionerTypeBase
  {
   public:
    MLPreconditioner(Teuchos::ParameterList& mllist);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override { return p_; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string get_parameter_list_name() const override { return "ML Parameters"; }

   private:
    Teuchos::ParameterList& mllist_;

    //! system of equations used for preconditioning used by P_ only
    Teuchos::RCP<Epetra_RowMatrix> pmatrix_;

    /// preconditioner
    Teuchos::RCP<Epetra_Operator> p_;
  };
}  // namespace CORE::LINEAR_SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
