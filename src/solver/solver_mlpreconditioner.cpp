/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class to ML preconditoner

\brief Declaration

\level 0

*/
/*----------------------------------------------------------------------*/

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"

#include "linalg_mlapi_operator.H"  // Michael's MLAPI based ML preconditioner
#include "linalg_utils_sparse_algebra_math.H"
#include "lib_dserror.H"

#include "solver_mlpreconditioner.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MLPreconditioner::MLPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist)
    : PreconditionerType(outfile), mllist_(mllist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MLPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
    if (A == NULL) dserror("CrsMatrix expected");

    // free old matrix first
    P_ = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    mllist_.remove("init smoother", false);

    // see whether we use standard ml or our own mlapi operator
    const bool domlapioperator = mllist_.get<bool>("LINALG::AMG_Operator", false);
    if (domlapioperator)
    {
      P_ = Teuchos::rcp(new LINALG::AMG_Operator(Pmatrix_, mllist_, true));
    }
    else
    {
      P_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Pmatrix_, mllist_, true));

      // for debugging ML
      // dynamic_cast<ML_Epetra::MultiLevelPreconditioner&>(*P_).PrintUnused(0);
    }
  }
}
