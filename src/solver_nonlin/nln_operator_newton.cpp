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
#include <Teuchos_TimeMonitor.hpp>

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
NLNSOL::NlnOperatorNewton::NlnOperatorNewton()
: linsolver_(Teuchos::null),
  linesearch_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNewton::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnOperatorNewton::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  if (FixedJacobian())
    jac_ = NlnProblem()->GetJacobianOperator();

  SetupLinearSolver();
  SetupLineSearch();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNewton::SetupLinearSolver()
{
  // get the solver number used for structural problems
  const int linsolvernumber = Params().get<int>("Newton: Linear Solver");

  // check if the solver ID is valid
  if (linsolvernumber == (-1))
    dserror("No valid linear solver defined!");

  linsolver_ = Teuchos::rcp(new LINALG::Solver(
      DRT::Problem::Instance()->SolverParams(linsolvernumber), Comm(),
      DRT::Problem::Instance()->ErrorFile()->Handle()));

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNewton::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ = linesearchfactory.Create(
      Params().sublist("Nonlinear Operator: Line Search"));

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNewton::ApplyInverse(const Epetra_MultiVector& f,
    Epetra_MultiVector& x) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnOperatorNewton::ApplyInverse");
  Teuchos::TimeMonitor monitor(*time);

  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // ---------------------------------------------------------------------------
  // initialize stuff for Newton loop
  // ---------------------------------------------------------------------------
  // solution increment vector
  Teuchos::RCP<Epetra_MultiVector> inc =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  // residual vector already including negative sign
  Teuchos::RCP<Epetra_MultiVector> rhs =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->ComputeF(x, *rhs);

  // some scalars
  int iter = 0; // iteration counter
  double steplength = 1.0; // line search parameter
  double fnorm2 = 1.0e+12; // residual L2 norm
  bool converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2); // convergence flag

  PrintIterSummary(iter, fnorm2);

  // ---------------------------------------------------------------------------
  // do a Newton-type iteration loop
  // ---------------------------------------------------------------------------
  while (ContinueIterations(iter, converged))
  {
    ++iter;

    // prepare linear solve
    NlnProblem()->ComputeJacobian();
    inc->PutScalar(0.0);

    // compute the Newton increment
    ComputeSearchDirection(rhs, inc, iter);

    // compute line search parameter
    steplength = ComputeStepLength(x, *rhs, *inc, fnorm2);

    // Iterative update
    err = x.Update(steplength, *inc, 1.0);
    if (err != 0) { dserror("Failed."); }

    // evaluate and check for convergence
    NlnProblem()->ComputeF(x, *rhs);
    converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2);

    PrintIterSummary(iter, fnorm2);
  }

  // return error code
  return (not CheckSuccessfulConvergence(iter, converged));
}

/*----------------------------------------------------------------------------*/
const int NLNSOL::NlnOperatorNewton::ComputeSearchDirection(
    Teuchos::RCP<Epetra_MultiVector>& rhs,
    Teuchos::RCP<Epetra_MultiVector>& inc, const int iter) const
{
  // error code for linear solver
  int linsolve_error = 0;

  // compute search direction with either fixed or most recent, updated jacobian
  if (FixedJacobian())
  {
    // do not refactor or reset since we want to reuse the LU decomposition
    linsolve_error = linsolver_->Solve(jac_, inc, rhs, false, false,
        Teuchos::null);
  }
  else
  {
    linsolve_error = linsolver_->Solve(NlnProblem()->GetJacobianOperator(), inc,
        rhs, true, iter == 1, Teuchos::null);
  }

  // test for failure of linear solver
  if (linsolve_error != 0) { dserror("Linear solver failed."); }

  // Is 'inc' a descent direction?
  double dotproduct = 0.0;
  inc->Dot(*rhs, &dotproduct);
  if (dotproduct <= 0.0) // ascent direction
    dserror("Search direction in %s is not a descent direction anymore.",
        Label());

  return linsolve_error;
}

/*----------------------------------------------------------------------------*/
const bool NLNSOL::NlnOperatorNewton::FixedJacobian() const
{
  return Params().get<bool>("Newton: Fixed Jacobian");
}

/*----------------------------------------------------------------------------*/
const double NLNSOL::NlnOperatorNewton::ComputeStepLength(
    const Epetra_MultiVector& x, const Epetra_MultiVector& f,
    const Epetra_MultiVector& inc, double fnorm2) const
{
  linesearch_->Init(NlnProblem(),
      Params().sublist("Nonlinear Operator: Line Search"), x, f, inc, fnorm2);
  linesearch_->Setup();
  return linesearch_->ComputeLSParam();
}
