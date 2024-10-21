// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linear_solver_preconditioner_krylovprojection.hpp"

#include "4C_linalg_projected_precond.hpp"

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::KrylovProjectionPreconditioner::KrylovProjectionPreconditioner(
    Teuchos::RCP<Core::LinearSolver::PreconditionerTypeBase> preconditioner,
    Teuchos::RCP<Core::LinAlg::KrylovProjector> projector)
    : preconditioner_(preconditioner), projector_(projector)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::KrylovProjectionPreconditioner::setup(bool create, Epetra_Operator* matrix,
    Core::LinAlg::MultiVector<double>* x, Core::LinAlg::MultiVector<double>* b)
{
  projector_->apply_pt(*b);

  // setup wrapped preconditioner
  preconditioner_->setup(create, matrix, x, b);

  // Wrap the linear operator of the contained preconditioner. This way the
  // actual preconditioner is called first and the projection is done
  // afterwards.

  p_ = Teuchos::make_rcp<Core::LinAlg::LinalgPrecondOperator>(
      preconditioner_->prec_operator(), true, projector_);
}

FOUR_C_NAMESPACE_CLOSE
