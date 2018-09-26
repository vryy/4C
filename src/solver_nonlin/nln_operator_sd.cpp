/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_sd.cpp

\brief Steepest-descent algorithm

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

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
#include "nln_operator_sd.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorSD::NlnOperatorSD() : linesearch_(Teuchos::null) { return; }

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorSD::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorSD::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  SetupLineSearch();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorSD::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ =
      linesearchfactory.Create(Configuration(), MyGetParameter<std::string>("line search"));

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorSD::ApplyInverse(const Epetra_MultiVector& f, Epetra_MultiVector& x) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorSD::ApplyInverse");
  Teuchos::TimeMonitor monitor(*time);

  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  // ---------------------------------------------------------------------------
  // initialize stuff for iteration loop
  // ---------------------------------------------------------------------------
  // solution increment vector
  Teuchos::RCP<Epetra_MultiVector> inc = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  // residual vector
  Teuchos::RCP<Epetra_MultiVector> rhs = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->ComputeF(x, *rhs);

  int iter = 0;                                                   // iteration counter
  double steplength = 1.0;                                        // line search parameter
  double fnorm2 = 1.0e+12;                                        // residual L2 norm
  bool converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2);  // convergence flag
  bool suffdecr = false;  // flag for sufficient decrease of line search

  // print initial state before iterating
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
    ComputeStepLength(x, *rhs, *inc, fnorm2, steplength, suffdecr);

    // Iterative update
    err = x.Update(steplength, *inc, 1.0);
    if (err != 0)
    {
      dserror("Failed.");
    }

    // compute current residual and check for convergence
    NlnProblem()->ComputeF(x, *rhs);
    converged = NlnProblem()->ConvergenceCheck(*rhs, fnorm2);

    PrintIterSummary(iter, fnorm2);
  }

  // ---------------------------------------------------------------------------
  // Finish ApplyInverse()
  // ---------------------------------------------------------------------------
  // determine error code
  NLNSOL::UTILS::OperatorStatus errorcode = ErrorCode(iter, converged, err);

  // write to output parameter list
  SetOutParameterIter(iter);
  SetOutParameterResidualNorm(fnorm2);
  SetOutParameterConverged(converged);
  SetOutParameterErrorCode(errorcode);

  // return error code
  return errorcode;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorSD::ComputeSearchDirection(
    const Epetra_MultiVector& rhs, Epetra_MultiVector& inc) const
{
  /* Since NlnProblem()->ComputeF() already returns the 'descending' residual
   * (= steepest descent direction), we just have to copy it into the increment
   * vector.
   */
  int err = inc.Update(1.0, rhs, 0.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  return err;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorSD::ComputeStepLength(const Epetra_MultiVector& x,
    const Epetra_MultiVector& f, const Epetra_MultiVector& inc, double fnorm2, double& lsparam,
    bool& suffdecr) const
{
  const std::string lslist = MyGetParameter<std::string>("line search");

  linesearch_->Init(NlnProblem(), Configuration(), lslist, x, f, inc, fnorm2);
  linesearch_->Setup();
  linesearch_->ComputeLSParam(lsparam, suffdecr);


  return;
}
