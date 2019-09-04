/*----------------------------------------------------------------------*/
/*! \file
\brief estimate condition number of a matrix

copied and adapted from
Trilinos/packages/ifpack/src/Ifpack_Condest.h and Ifpack_Condest.cpp

\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/


#include "AztecOO.h"
#include "AztecOO_StatusTestResNorm.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"

#include "linalg_solver.H"
#include "linalg_sparsematrix.H"
#include <Teuchos_TimeMonitor.hpp>



double LINALG::Condest(
    SparseMatrix& Matrix, const Ifpack_CondestType CT, const int MaxIters, const double Tol)
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Condest");

  double ConditionNumberEstimate = -1.0;

  if (CT == Ifpack_Cheap)
  {
    // Create a vector with all values equal to one
    Epetra_Vector Ones(Matrix.OperatorDomainMap());
    Ones.PutScalar(1.0);
    // Create the vector of results
    Epetra_Vector OnesResult(Matrix.OperatorRangeMap());
    // Compute the effect of the solve on the vector of ones
    IFPACK_CHK_ERR(Matrix.ApplyInverse(Ones, OnesResult));
    // Make all values non-negative
    IFPACK_CHK_ERR(OnesResult.Abs(OnesResult));
    // Get the maximum value across all processors
    IFPACK_CHK_ERR(OnesResult.MaxValue(&ConditionNumberEstimate));
  }
  else if (CT == Ifpack_CG)
  {
    Epetra_Vector LHS(Matrix.OperatorDomainMap());
    LHS.PutScalar(0.0);
    Epetra_Vector RHS(Matrix.OperatorRangeMap());
    RHS.Random();
    Epetra_LinearProblem Problem;
    Problem.SetOperator(Matrix.EpetraMatrix().get());
    Problem.SetLHS(&LHS);
    Problem.SetRHS(&RHS);

    AztecOO Solver(Problem);
    Solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
    Solver.SetAztecOption(AZ_output, AZ_none);
    // Solver.SetAztecOption(AZ_output,10);
    Solver.Iterate(MaxIters, Tol);

    const double* status = Solver.GetAztecStatus();
    ConditionNumberEstimate = status[AZ_condnum];
  }
  else if (CT == Ifpack_GMRES)
  {
    Epetra_Vector LHS(Matrix.OperatorDomainMap());
    LHS.PutScalar(0.0);
    Epetra_Vector RHS(Matrix.OperatorRangeMap());
    RHS.Random();
    Epetra_LinearProblem Problem;
    Problem.SetOperator(Matrix.EpetraMatrix().get());
    Problem.SetLHS(&LHS);
    Problem.SetRHS(&RHS);

    AztecOO Solver(Problem);
    Solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
    Solver.SetAztecOption(AZ_output, AZ_none);
    // Solver.SetAztecOption(AZ_output,10);

    // the following can be problematic for large problems,
    // but any restart would destroy useful information about
    // the condition number.
    const int iterlimit = 150;
    int MaxIters_mod = MaxIters;
    if (MaxIters > iterlimit and Matrix.Comm().MyPID() == 0)
    {
      //      MaxIters_mod = iterlimit;
      std::cout << std::endl << "Krylov space size: " << MaxIters_mod << std::endl;
      std::cout << "Warning! Large MaxIters means a large Krylov subspace for GMRES -> you might "
                   "run out of memory."
                << std::endl;
      std::cout << "Continuing with MaxIters = " << MaxIters_mod << std::endl;
    }
    Solver.SetAztecOption(AZ_kspace, MaxIters_mod);  // Krylov space is set to iteration number !!!
    Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
    Solver.Iterate(MaxIters_mod, Tol);

    const double* status = Solver.GetAztecStatus();
    ConditionNumberEstimate = status[AZ_condnum];
  }

  return (ConditionNumberEstimate);
}
