/*!----------------------------------------------------------------------

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
*----------------------------------------------------------------------*/

#include "../linalg/linalg_krylov_projector.H"
#include "../linalg/linalg_projected_operator.H"  // for LINALG::LinalgProjectedOperator
#include "../linalg/linalg_projected_precond.H"   // for LINALG::LinalgPrecondOperator

#include "solver_krylovprojectionpreconditioner.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::KrylovProjectionPreconditioner::KrylovProjectionPreconditioner(FILE* outfile,
    Teuchos::RCP<LINALG::SOLVER::PreconditionerType> preconditioner,
    Teuchos::RCP<LINALG::KrylovProjector> projector)
    : LINALG::SOLVER::PreconditionerType(outfile),
      preconditioner_(preconditioner),
      projector_(projector)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovProjectionPreconditioner::Setup(
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
void LINALG::SOLVER::KrylovProjectionPreconditioner::Finish(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  preconditioner_->Finish(matrix, x, b);
}
