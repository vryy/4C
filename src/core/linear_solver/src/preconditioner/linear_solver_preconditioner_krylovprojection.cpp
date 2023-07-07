/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

 \level 1

*/
/*----------------------------------------------------------------------*/

#include "linalg_projected_precond.H"

#include "linear_solver_preconditioner_krylovprojection.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner::KrylovProjectionPreconditioner(FILE* outfile,
    Teuchos::RCP<CORE::LINEAR_SOLVER::PreconditionerType> preconditioner,
    Teuchos::RCP<LINALG::KrylovProjector> projector)
    : CORE::LINEAR_SOLVER::PreconditionerType(outfile),
      preconditioner_(preconditioner),
      projector_(projector)
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

  A_ = Teuchos::rcp(new LINALG::LinalgProjectedOperator(Teuchos::rcp(A, false), true, projector_));

  P_ = Teuchos::rcp(new LINALG::LinalgPrecondOperator(
      Teuchos::rcp(preconditioner_->PrecOperator(), false), true, projector_));

  SetupLinearProblem(&*A_, lp.GetLHS(), lp.GetRHS());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner::Finish(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  preconditioner_->Finish(matrix, x, b);
}
