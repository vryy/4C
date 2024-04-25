/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

 \level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_linear_solver_preconditioner_krylovprojection.hpp"

#include "4C_linalg_projected_precond.hpp"

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner::KrylovProjectionPreconditioner(
    Teuchos::RCP<CORE::LINEAR_SOLVER::PreconditionerTypeBase> preconditioner,
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
    : preconditioner_(preconditioner), projector_(projector)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  projector_->ApplyPT(*b);

  // setup wrapped preconditioner
  preconditioner_->Setup(create, matrix, x, b);

  // Wrap the linear operator of the contained preconditioner. This way the
  // actual preconditioner is called first and the projection is done
  // afterwards.

  Epetra_LinearProblem& lp = preconditioner_->LinearProblem();
  Epetra_Operator* A = lp.GetOperator();

  a_ = Teuchos::rcp(
      new CORE::LINALG::LinalgProjectedOperator(Teuchos::rcp(A, false), true, projector_));

  p_ = Teuchos::rcp(
      new CORE::LINALG::LinalgPrecondOperator(preconditioner_->PrecOperator(), true, projector_));

  SetupLinearProblem(&*a_, lp.GetLHS(), lp.GetRHS());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner::Finish(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  preconditioner_->Finish(matrix, x, b);
}

FOUR_C_NAMESPACE_CLOSE
