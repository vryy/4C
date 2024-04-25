/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 1


*----------------------------------------------------------------------*/

#include "4C_linalg_projected_operator.hpp"

#include "4C_linalg_krylov_projector.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_MultiVector.h>

FOUR_C_NAMESPACE_OPEN

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
CORE::LINALG::LinalgProjectedOperator::LinalgProjectedOperator(Teuchos::RCP<Epetra_Operator> A,
    bool project, Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
    : project_(project), a_(A), projector_(projector)
{
  if (project_ && (projector == Teuchos::null))
    FOUR_C_THROW("Kernel projection enabled but got no projector object");

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
  ierr = a_->Apply(X, Y);

  // if necessary, project out matrix kernel
  if (project_)
  {
    // int ierr2=0;
    /*ierr2 = */
    projector_->ApplyPT(Y);
  }

  return (ierr);
}  // CORE::LINALG::LinalgProjectedOperator::Apply

FOUR_C_NAMESPACE_CLOSE
