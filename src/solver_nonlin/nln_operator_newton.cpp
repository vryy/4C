/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_newton.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
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
#include "nln_operator_newton.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperatorNewton::NlnOperatorNewton()
: linsolver_(Teuchos::null),
  linesearch_(Teuchos::null),
  maxiter_(1),
  fixedjacobian_(false)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of the algorithm  / operator */
void NLNSOL::NlnOperatorNewton::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // ---------------------------------------------------------------------------
  // initialize member variables from parameter list
  // ---------------------------------------------------------------------------
  if (Params().isParameter("Newton: Max Iter"))
    maxiter_ = Params().get<int>("Newton: Max Iter");

  if (Params().isParameter("Newton: Fixed Jacobian"))
    fixedjacobian_ = Params().get<bool>("Newton: Fixed Jacobian");

  if (fixedjacobian_)
    jac_ = NlnProblem()->GetJacobianOperator();

  // ---------------------------------------------------------------------------

  SetupLinearSolver();
  SetupLineSearch();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of linear solver */
void NLNSOL::NlnOperatorNewton::SetupLinearSolver()
{
  // get the solver number used for structural problems
  const int linsolvernumber = Params().get<int>("Newton: Linear Solver");

  // check if the solver ID is valid
  if (linsolvernumber == (-1))
    dserror("No valid linear solver defined!");

  linsolver_ = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                               Comm(),
                                               DRT::Problem::Instance()->ErrorFile()->Handle()));

  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of line search */
void NLNSOL::NlnOperatorNewton::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ = linesearchfactory.Create(Params().sublist("Newton: Line Search"));

  return;
}

/*----------------------------------------------------------------------------*/
/* Apply the preconditioner */
int NLNSOL::NlnOperatorNewton::ApplyInverse(const Epetra_MultiVector& f,
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
  double steplength = 1.0; // line search parameter
  double fnorm2 = 1.0e+12; // residual L2 norm
  bool converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2); // convergence flag

  if (Comm().MyPID() == 0 and Params().get<bool>("Newton: Print Iterations"))
    IO::cout << "Newton-type iteration " << iter << ": res-norm = " << fnorm2 << IO::endl;

  // ---------------------------------------------------------------------------
  // do a Newton-type iteration loop
  // ---------------------------------------------------------------------------
  while (iter < GetMaxIter() && not converged)
  {
    ++iter;

    // prepare linear solve
    rhs->Scale(-1.0);
    inc->PutScalar(0.0);

    // compute the Newton increment
    ComputeSearchDirection(inc, rhs, iter);

    // line search
    linesearch_->Init(NlnProblem(), Params().sublist("Newton: Line Search"), x, *inc, fnorm2);
    linesearch_->Setup();
    steplength = linesearch_->ComputeLSParam();

    // Iterative update
    err = x.Update(steplength, *inc, 1.0);
    if (err != 0) { dserror("Failed."); }

    // evaluate and check for convergence
    NlnProblem()->Evaluate(x, *rhs);

    converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2);

    if (Params().get<bool>("Newton: Print Iterations") and Comm().MyPID() == 0)
      IO::cout << "Newton-type iteration " << iter
               << ": res-norm = " << std::setprecision(6) << fnorm2
               << IO::endl;
  }

  // ---------------------------------------------------------------------------
  // check successful convergence
  // ---------------------------------------------------------------------------
  if (iter <= GetMaxIter() && converged)
    return 0;
  else
  {
    if (IsSolver())
      dserror("Newton-type algorithm did not converge in %d iterations.", iter);
  }

  // suppose non-convergence
  return 1;
}

/*----------------------------------------------------------------------------*/
/* Compute the search direction */
const int NLNSOL::NlnOperatorNewton::ComputeSearchDirection(
    Teuchos::RCP<Epetra_MultiVector>& inc,
    Teuchos::RCP<Epetra_MultiVector>& rhs,
    const int iter) const
{
  // error code for linear solver
  int linsolve_error = 0;

  // compute search direction with either fixed or most recent, updated jacobian
  if (fixedjacobian_)
  {
    linsolve_error = linsolver_->Solve(jac_, inc, rhs, true, iter==1, Teuchos::null);
  }
  else
  {
    linsolve_error = linsolver_->Solve(NlnProblem()->GetJacobianOperator(), inc, rhs, true, iter==1, Teuchos::null);
  }

  // test for failure of linear solver
  if (linsolve_error != 0) { dserror("Linear solver failed."); }

  return linsolve_error;
}
