/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class to ML preconditoner

\brief Declaration

\level 0

*/
/*----------------------------------------------------------------------*/

#include "4C_linear_solver_preconditioner_ml.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"

#include <ml_common.h>
#include <ml_epetra.h>
#include <ml_epetra_operator.h>
#include <ml_epetra_utils.h>
#include <ml_include.h>
#include <ml_MultiLevelPreconditioner.h>

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::MLPreconditioner::MLPreconditioner(Teuchos::ParameterList& mllist)
    : mllist_(mllist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::MLPreconditioner::setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  if (create)
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
    if (A == nullptr) FOUR_C_THROW("CrsMatrix expected");

    // free old matrix first
    p_ = Teuchos::null;
    pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    mllist_.remove("init smoother", false);

    p_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*pmatrix_, mllist_, true));
  }
}

FOUR_C_NAMESPACE_CLOSE
