/*!----------------------------------------------------------------------

\brief Implementation

\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

*----------------------------------------------------------------------*/

#include "linalg_projected_precond.H"
#include "Epetra_MultiVector.h"
#include "linalg_krylov_projector.H"
#include "../drt_lib/drt_dserror.H"

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
LINALG::LinalgPrecondOperator::LinalgPrecondOperator(Teuchos::RCP<Epetra_Operator> precond,
    bool project, Teuchos::RCP<LINALG::KrylovProjector> projector)
    : project_(project), precond_(precond), projector_(projector)
{
  if (project_ && (projector == Teuchos::null))
    dserror("Kernel projection enabled but got no projector object");

  return;
}  // LINALG::LinalgPrecondOperator::LinalgPrecondOperator

/* --------------------------------------------------------------------
                          Destructor
   -------------------------------------------------------------------- */
LINALG::LinalgPrecondOperator::~LinalgPrecondOperator()
{
  return;
}  // LINALG::LinalgPrecondOperator::~KrylovProjector

/* --------------------------------------------------------------------
                    (Modified) ApplyInverse call
   -------------------------------------------------------------------- */
int LINALG::LinalgPrecondOperator::ApplyInverse(
    const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  int ierr = 0;
  // Apply the inverse preconditioner to get new basis vector for the
  // Krylov space
  ierr = precond_->ApplyInverse(X, Y);

  // if necessary, project out matrix kernel to maintain well-posedness
  // of problem
  if (project_)
  {
    projector_->ApplyP(Y);
  }

  return (ierr);
}  // LINALG::LinalgPrecondOperator::ApplyInverse
