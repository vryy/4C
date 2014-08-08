/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_sd.cpp

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
#include "nln_operator_sd.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperatorSD::NlnOperatorSD()
: linesearch_(Teuchos::null),
  maxiter_(1)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of the algorithm  / operator */
void NLNSOL::NlnOperatorSD::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // ---------------------------------------------------------------------------
  // initialize member variables from parameter list
  // ---------------------------------------------------------------------------
  if (Params().isParameter("SD: Max Iter"))
    maxiter_ = Params().get<int>("SD: Max Iter");

  // ---------------------------------------------------------------------------

  SetupLineSearch();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of line search */
void NLNSOL::NlnOperatorSD::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ = linesearchfactory.Create(Params().sublist("SD: Line Search"));

  return;
}

/*----------------------------------------------------------------------------*/
/* Apply the preconditioner */
int NLNSOL::NlnOperatorSD::ApplyInverse(const Epetra_MultiVector& f,
    Epetra_MultiVector& x) const
{
  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // ---------------------------------------------------------------------------
  // initialize stuff for iteration loop
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

  if (Params().get<bool>("SD: Print Iterations"))
    PrintIterSummary(iter, fnorm2);

  // ---------------------------------------------------------------------------
  // iteration loop
  // ---------------------------------------------------------------------------
  while (iter < GetMaxIter() && not converged)
  {
    ++iter;

    // compute the search direction
    ComputeSearchDirection(*rhs, *inc);

    // line search
    linesearch_->Init(NlnProblem(), Params().sublist("SD: Line Search"), x, *inc, fnorm2);
    linesearch_->Setup();
    steplength = linesearch_->ComputeLSParam();

    // Iterative update
    err = x.Update(steplength, *inc, 1.0);
    if (err != 0) { dserror("Failed."); }

    // evaluate and check for convergence
    NlnProblem()->Evaluate(x, *rhs);

    converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2);

    if (Params().get<bool>("SD: Print Iterations"))
      PrintIterSummary(iter, fnorm2);
  }

  // ---------------------------------------------------------------------------
  // check successful convergence
  // ---------------------------------------------------------------------------
  if (iter <= GetMaxIter() && converged)
    return 0;
  else
  {
    if (IsSolver())
      dserror("Steepest descent algorithm did not converge in %d iterations.", iter);
  }

  // suppose non-convergence
  return 1;
}

/*----------------------------------------------------------------------------*/
/* Compute the search direction */
const int NLNSOL::NlnOperatorSD::ComputeSearchDirection(
    const Epetra_MultiVector& rhs,
    Epetra_MultiVector& inc
    ) const
{
  // search direction = negative residual (gradient)
  int err = inc.Update(-1.0, rhs, 0.0);
  if (err != 0) { dserror("Update failed."); }

  return err;
}
