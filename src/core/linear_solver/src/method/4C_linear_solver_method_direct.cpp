// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linear_solver_method_direct.hpp"

#include "4C_linalg_krylov_projector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Amesos_Klu.h>
#include <Amesos_Superludist.h>
#include <Amesos_Umfpack.h>
#include <Epetra_LinearProblem.h>

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
Core::LinearSolver::DirectSolver<MatrixType, VectorType>::DirectSolver(std::string solvertype)
    : solvertype_(solvertype),
      factored_(false),
      solver_(nullptr),
      reindexer_(nullptr),
      projector_(nullptr)
{
  linear_problem_ = std::make_shared<Epetra_LinearProblem>();
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void Core::LinearSolver::DirectSolver<MatrixType, VectorType>::setup(
    std::shared_ptr<MatrixType> matrix, std::shared_ptr<VectorType> x,
    std::shared_ptr<VectorType> b, const bool refactor, const bool reset,
    std::shared_ptr<Core::LinAlg::KrylovProjector> projector)
{
  std::shared_ptr<Epetra_CrsMatrix> crsA = std::dynamic_pointer_cast<Epetra_CrsMatrix>(matrix);

  // 1. merge the block system matrix into a standard sparse matrix if necessary
  if (!crsA)
  {
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> Ablock =
        std::dynamic_pointer_cast<Core::LinAlg::BlockSparseMatrixBase>(matrix);

    int matrixDim = Ablock->full_range_map().NumGlobalElements();
    if (matrixDim > 50000)
      std::cout << "\n WARNING: Direct linear solver is merging matrix, this is very expensive! \n";

    std::shared_ptr<Core::LinAlg::SparseMatrix> Ablock_merged = Ablock->merge();
    crsA = Ablock_merged->epetra_matrix();
  }

  // 2. project the linear system if close to being singular and set the final matrix and vectors
  projector_ = projector;
  if (projector_ != nullptr)
  {
    Core::LinAlg::SparseMatrix A_view(crsA, Core::LinAlg::View);
    std::shared_ptr<Core::LinAlg::SparseMatrix> A2 = projector_->project(A_view);

    crsA = A2->epetra_matrix();
    projector_->apply_pt(*b);
  }

  x_ = x;
  b_ = b;
  a_ = crsA;

  // 3. Do a GID reindexing of the overall problem and create the direct solver
  linear_problem_->SetRHS(b_->get_ptr_of_Epetra_MultiVector().get());
  linear_problem_->SetLHS(x_->get_ptr_of_Epetra_MultiVector().get());
  linear_problem_->SetOperator(a_.get());

  if (reindexer_ and not(reset or refactor)) reindexer_->fwd();

  if (reset or refactor or not is_factored())
  {
    reindexer_ = std::make_shared<EpetraExt::LinearProblem_Reindex2>(nullptr);

    if (solvertype_ == "umfpack")
    {
      solver_ = std::make_shared<Amesos_Umfpack>((*reindexer_)(*linear_problem_));
    }
    else if (solvertype_ == "superlu")
    {
      solver_ = std::make_shared<Amesos_Superludist>((*reindexer_)(*linear_problem_));
    }
    else
    {
      solver_ = std::make_shared<Amesos_Klu>((*reindexer_)(*linear_problem_));
    }

    factored_ = false;
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
int Core::LinearSolver::DirectSolver<MatrixType, VectorType>::solve()
{
  if (not is_factored())
  {
    solver_->SymbolicFactorization();
    solver_->NumericFactorization();
    factored_ = true;
  }

  solver_->Solve();

  if (projector_ != nullptr) projector_->apply_p(*x_);

  return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// explicit initialization
template class Core::LinearSolver::DirectSolver<Epetra_Operator, Core::LinAlg::MultiVector<double>>;

FOUR_C_NAMESPACE_CLOSE
