/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 1

*----------------------------------------------------------------------*/

#include "4C_linalg_projected_precond.hpp"

#include "4C_linalg_krylov_projector.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_MultiVector.h>

FOUR_C_NAMESPACE_OPEN

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
Core::LinAlg::LinalgPrecondOperator::LinalgPrecondOperator(Teuchos::RCP<Epetra_Operator> precond,
    bool project, Teuchos::RCP<Core::LinAlg::KrylovProjector> projector)
    : project_(project), precond_(precond), projector_(projector)
{
  if (project_ && (projector == Teuchos::null))
    FOUR_C_THROW("Kernel projection enabled but got no projector object");

  return;
}  // Core::LinAlg::LinalgPrecondOperator::LinalgPrecondOperator


/* --------------------------------------------------------------------
                    (Modified) ApplyInverse call
   -------------------------------------------------------------------- */
int Core::LinAlg::LinalgPrecondOperator::ApplyInverse(
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
}  // Core::LinAlg::LinalgPrecondOperator::ApplyInverse

FOUR_C_NAMESPACE_CLOSE
