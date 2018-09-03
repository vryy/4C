/*!----------------------------------------------------------------------
\file  linalg_projected_operator.cpp

\brief Implementation

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

*----------------------------------------------------------------------*/

#include "linalg_projected_operator.H"
#include "Epetra_MultiVector.h"
#include "linalg_krylov_projector.H"
#include "../drt_lib/drt_dserror.H"

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
LINALG::LinalgProjectedOperator::LinalgProjectedOperator(
    Teuchos::RCP<Epetra_Operator> A, bool project, Teuchos::RCP<LINALG::KrylovProjector> projector)
    : project_(project), A_(A), projector_(projector)
{
  if (project_ && (projector == Teuchos::null))
    dserror("Kernel projection enabled but got no projector object");

  return;
}  // LINALG::LinalgProjectedOperator::LinalgProjectedOperator

/* --------------------------------------------------------------------
                           Destructor
   -------------------------------------------------------------------- */
LINALG::LinalgProjectedOperator::~LinalgProjectedOperator()
{
  return;
}  // LINALG::LinalgProjectedOperator::~LinalgProjectedOperator

/* --------------------------------------------------------------------
                      (Modified) Apply call
   -------------------------------------------------------------------- */
int LINALG::LinalgProjectedOperator::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
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
}  // LINALG::LinalgProjectedOperator::Apply
