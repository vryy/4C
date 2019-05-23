/*!----------------------------------------------------------------------

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
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
