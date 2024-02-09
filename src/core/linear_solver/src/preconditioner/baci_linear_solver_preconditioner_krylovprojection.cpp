/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

 \level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_linear_solver_preconditioner_krylovprojection.hpp"

#include "baci_linalg_projected_precond.hpp"

BACI_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner::KrylovProjectionPreconditioner(
    Teuchos::RCP<CORE::LINEAR_SOLVER::PreconditionerType> preconditioner,
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

  A_ = Teuchos::rcp(
      new CORE::LINALG::LinalgProjectedOperator(Teuchos::rcp(A, false), true, projector_));

  P_ = Teuchos::rcp(
      new CORE::LINALG::LinalgPrecondOperator(preconditioner_->PrecOperator(), true, projector_));

  SetupLinearProblem(&*A_, lp.GetLHS(), lp.GetRHS());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner::Finish(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  preconditioner_->Finish(matrix, x, b);
}

BACI_NAMESPACE_CLOSE
