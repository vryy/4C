/*----------------------------------------------------------------------*/
/*!
\file nln_operator_quasinewton.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>

// standard
#include <iostream>
#include <vector>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "linesearch_base.H"
#include "linesearch_factory.H"
#include "nln_operator_quasinewton.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperatorQuasiNewton::NlnOperatorQuasiNewton()
: linsolver_(Teuchos::null),
  linesearch_(Teuchos::null),
  maxiter_(1)
{
  return;
}

/*----------------------------------------------------------------------*/
/* Setup of the algorithm  / operator */
void NLNSOL::NlnOperatorQuasiNewton::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  maxiter_ = Params().get<int>("Quasi Newton: Max Iter");

  SetupLinearSolver();
  SetupLineSearch();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------*/
/* Setup of linear solver */
void NLNSOL::NlnOperatorQuasiNewton::SetupLinearSolver()
{
  // get the solver number used for structural problems
  const int linsolvernumber = Params().get<int>("Quasi Newton: Linear Solver");

  // check if the solver ID is valid
  if (linsolvernumber == (-1))
    dserror("No valid linear solver defined!");

  linsolver_ = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                               Comm(),
                                               DRT::Problem::Instance()->ErrorFile()->Handle()));

  return;
}

/*----------------------------------------------------------------------*/
/* Setup of line search */
void NLNSOL::NlnOperatorQuasiNewton::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ = linesearchfactory.Create(Params().sublist("Quasi Newton: Line Search"));

  return;
}

/*----------------------------------------------------------------------*/
/* Apply the preconditioner */
int NLNSOL::NlnOperatorQuasiNewton::ApplyInverse(const Epetra_MultiVector& f,
    Epetra_MultiVector& x) const
{
  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // ---------------------------------------------------------------------------
  // initialize stuff for Newton loop
  // ---------------------------------------------------------------------------
  // solution increment vector
  Teuchos::RCP<Epetra_MultiVector> inc = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  // residual vector
  Teuchos::RCP<Epetra_MultiVector> rhs = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->Evaluate(x, *rhs);

  // some scalars
  int iter = 0; // iteration counter
  int linsolve_error = 0; // error code for linear solver
  double steplength = 1.0; // line search parameter
  double fnorm2 = 1.0e+12; // residual L2 norm
  bool converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2); // convergence flag

  if (Comm().MyPID() == 0 and Params().get<bool>("Quasi Newton: Print Iterations"))
    IO::cout << "Quasi-Newton iteration " << iter << ": res-norm = " << fnorm2 << IO::endl;

  // ---------------------------------------------------------------------------
  // do a Quasi-Newton scheme with fixed Jacobian
  // ---------------------------------------------------------------------------
  while (iter < GetMaxIter() && not converged)
  {
    ++iter;

    // prepare linear solve
    rhs->Scale(-1.0);
    inc->PutScalar(0.0);

    // do the linear solve (no re-factorization since matrix has not changed)
    linsolve_error = linsolver_->Solve(NlnProblem()->GetJacobianOperator(), inc, rhs, iter==1, iter==1, Teuchos::null);
    if (linsolve_error != 0) { dserror("Linear solver failed."); }

    // line search
    linesearch_->Init(NlnProblem(), Params().sublist("Quasi Newton: Line Search"), x, *inc, fnorm2);
    linesearch_->Setup();
    steplength = linesearch_->ComputeLSParam();

    // Iterative update
    err = x.Update(steplength, *inc, 1.0);
    if (err != 0) { dserror("Failed."); }

    // evaluate and check for convergence
    NlnProblem()->Evaluate(x, *rhs);
    converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2);

    if (Params().get<bool>("Quasi Newton: Print Iterations") and Comm().MyPID() == 0)
      IO::cout << "Quasi-Newton iteration " << iter << ": res-norm = " << fnorm2 << IO::endl;
  }

  // ---------------------------------------------------------------------------
  // check successful convergence
  // ---------------------------------------------------------------------------
  if (iter < GetMaxIter() && converged)
    return 0;
  else
  {
    if (IsSolver())
      dserror("Quasi-Newton did not converge in %d iterations.", iter);
  }

  // suppose non-convergence
  return 1;
}
