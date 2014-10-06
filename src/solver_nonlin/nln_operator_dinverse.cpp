/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_dinverse.cpp

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
#include <Epetra_CrsMatrix.h>
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
#include "nln_operator_dinverse.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorDInverse::NlnOperatorDInverse()
: linesearch_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorDInverse::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  SetupLineSearch();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorDInverse::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ = linesearchfactory.Create(Params().sublist("Inverse Diagonal: Line Search"));

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorDInverse::ApplyInverse(const Epetra_MultiVector& f,
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
  Teuchos::RCP<Epetra_MultiVector> inc =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  // residual vector
  Teuchos::RCP<Epetra_MultiVector> rhs =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->ComputeF(x, *rhs);
  NlnProblem()->ComputeJacobian();

  // some scalars
  int iter = 0; // iteration counter
  double steplength = 1.0; // line search parameter
  double fnorm2 = 1.0e+12; // residual L2 norm
  bool converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2); // convergence flag

  PrintIterSummary(iter, fnorm2);

  // ---------------------------------------------------------------------------
  // iteration loop
  // ---------------------------------------------------------------------------
  while (ContinueIterations(iter, converged))
  {
    ++iter;

    // compute the search direction
    ComputeSearchDirection(*rhs, *inc);

    // line search
    steplength = ComputeStepLength(x, *inc, fnorm2);

    // Iterative update
    err = x.Update(steplength, *inc, 1.0);
    if (err != 0) { dserror("Failed."); }

    // compute current residual and check for convergence
    NlnProblem()->ComputeF(x, *rhs);
    NlnProblem()->ComputeJacobian();
    converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2);

    PrintIterSummary(iter, fnorm2);
  }

  // return error code
  return (not CheckSuccessfulConvergence(iter, converged));
}

/*----------------------------------------------------------------------------*/
const int NLNSOL::NlnOperatorDInverse::ComputeSearchDirection(
    const Epetra_MultiVector& rhs, Epetra_MultiVector& inc) const
{
  int err = 0;

  Teuchos::RCP<Epetra_Vector> stiffdiagvec = Teuchos::rcp(
      new Epetra_Vector(NlnProblem()->GetJacobianOperator()->OperatorRangeMap(),
          true));
  Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
      NlnProblem()->GetJacobianOperator(), true)->ExtractDiagonalCopy(
      *stiffdiagvec);
  Teuchos::RCP<LINALG::SparseMatrix> stiffdiag = LINALG::Eye(
      NlnProblem()->GetJacobianOperator()->OperatorRangeMap());
  err = stiffdiag->ReplaceDiagonalValues(*stiffdiagvec);
  if (err != 0) { dserror("RepalceDiagonalValues failed."); }

  err = stiffdiag->ApplyInverse(rhs, inc);
  if (err != 0) { dserror("Update failed."); }

  inc.Scale(-1.0);

  return err;
}

/*----------------------------------------------------------------------------*/
const double NLNSOL::NlnOperatorDInverse::ComputeStepLength(
    const Epetra_MultiVector& x, const Epetra_MultiVector& inc,
    double fnorm2) const
{
  linesearch_->Init(NlnProblem(), Params().sublist("Inverse Diagonal: Line Search"), x, inc,
      fnorm2);
  linesearch_->Setup();
  return linesearch_->ComputeLSParam();
}
