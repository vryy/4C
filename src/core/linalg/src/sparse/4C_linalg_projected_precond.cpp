// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_projected_precond.hpp"

#include "4C_linalg_krylov_projector.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_exceptions.hpp"


FOUR_C_NAMESPACE_OPEN

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
Core::LinAlg::LinalgPrecondOperator::LinalgPrecondOperator(std::shared_ptr<Epetra_Operator> precond,
    bool project, std::shared_ptr<Core::LinAlg::KrylovProjector> projector)
    : project_(project), precond_(precond), projector_(projector)
{
  if (project_ && (projector == nullptr))
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
    VectorView Y_view(Y);
    projector_->apply_p(Y_view);
  }

  return (ierr);
}  // Core::LinAlg::LinalgPrecondOperator::ApplyInverse

FOUR_C_NAMESPACE_CLOSE
