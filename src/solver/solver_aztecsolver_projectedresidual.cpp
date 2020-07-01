/*----------------------------------------------------------------------*/
/*! \file


\brief Declaration
\level 1

Created on: Mar 4, 2013

*---------------------------------------------------------------------------*/

#include "solver_aztecsolver_projectedresidual.H"
#include "../linalg/linalg_krylov_projector.H"

/* ====================================================================
    public
   ==================================================================== */

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
AztecOO_StatusTestProjResNorm::AztecOO_StatusTestProjResNorm(const Epetra_Operator& Operator,
    const Epetra_Vector& LHS, const Epetra_Vector& RHS,
    const Teuchos::RCP<LINALG::KrylovProjector>& projector, double Tolerance)
    : AztecOO_StatusTestResNorm(Operator, LHS, RHS, Tolerance)
{
  projector_ = projector;
}

/* --------------------------------------------------------------------
             Convergence check with residual projection
   -------------------------------------------------------------------- */
AztecOO_StatusType AztecOO_StatusTestProjResNorm::CheckStatus(int CurrentIter,
    Epetra_MultiVector* CurrentResVector, double CurrentResNormEst, bool SolutionUpdated)
{
  // project residual - the only difference to original method
  projector_->ApplyPT(*CurrentResVector);

  return AztecOO_StatusTestResNorm::CheckStatus(
      CurrentIter, CurrentResVector, -1.0, SolutionUpdated);
}
