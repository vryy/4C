/*
 * solver_preconditionertype.cpp
 *
 *  Created on: Jul 4, 2011
 *      Author: wiesner
 */

#include "solver_preconditionertype.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::PreconditionerType::PreconditionerType( FILE * outfile )
  : outfile_( outfile )
{
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::PreconditionerType::SetupLinearProblem( Epetra_Operator * matrix,
                                                             Epetra_MultiVector * x,
                                                             Epetra_MultiVector * b )
{
  lp_.SetOperator(matrix);
  lp_.SetLHS(x);
  lp_.SetRHS(b);
}
