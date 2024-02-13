/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 1


*----------------------------------------------------------------------*/

#include "baci_linalg_projected_operator.hpp"

#include "baci_linalg_krylov_projector.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_MultiVector.h>

BACI_NAMESPACE_OPEN

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
CORE::LINALG::LinalgProjectedOperator::LinalgProjectedOperator(Teuchos::RCP<Epetra_Operator> A,
    bool project, Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
    : project_(project), A_(A), projector_(projector)
{
  if (project_ && (projector == Teuchos::null))
    dserror("Kernel projection enabled but got no projector object");

  return;
}  // CORE::LINALG::LinalgProjectedOperator::LinalgProjectedOperator


/* --------------------------------------------------------------------
                      (Modified) Apply call
   -------------------------------------------------------------------- */
int CORE::LINALG::LinalgProjectedOperator::Apply(
    const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  int ierr = 0;

  // Apply the operator
  ierr = A_->Apply(X, Y);

  // if necessary, project out matrix kernel
  if (project_)
  {
    // int ierr2=0;
    /*ierr2 = */
    projector_->ApplyPT(Y);
  }

  return (ierr);
}  // CORE::LINALG::LinalgProjectedOperator::Apply

BACI_NAMESPACE_CLOSE
