/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 1
Created on: Jul 4, 2011
*----------------------------------------------------------------------*/

#include "solver_preconditionertype.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::PreconditionerType::PreconditionerType(FILE* outfile) : outfile_(outfile) {}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::PreconditionerType::SetupLinearProblem(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  lp_.SetOperator(matrix);
  lp_.SetLHS(x);
  lp_.SetRHS(b);
}
