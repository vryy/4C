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
#include "nln_utils.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorNewton::NlnOperatorNewton()
    : jacevery_(1), linsolver_(Teuchos::null), linesearch_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNewton::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorNewton::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  {
    jacevery_ = MyGetParameter<int>("Newton: update jacobian every");

    switch (jacevery_)
    {
      case 0:
      {
        if (getVerbLevel() > Teuchos::VERB_NONE)
        {
          *getOStream() << LabelShort() << " with fixed Jacobian: Chord's method" << std::endl;
        }
        break;
      }
      case 1:
      {
        if (getVerbLevel() > Teuchos::VERB_NONE)
        {
          *getOStream() << LabelShort() << " with updated Jacobian: Full Newton" << std::endl;
        }
        break;
      }
      default:
      {
        dserror(
            "Updating the Jacobian every %i iterations is not allowed. "
            "Choose '0' for Chord's method and '1' for Full Newton.",
            jacevery_);
        break;
      }
    }
  }

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
  const int linsolvernumber = MyGetParameter<int>("Newton: Linear Solver");

  switch (linsolvernumber)
  {
    case 0:  // use solver passed in from Baci
    {
      linsolver_ = BaciLinearSolver();
      break;
    }
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    {
      linsolver_ =
          Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
              Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
      break;
    }
    default:
    {
      dserror(
          "No valid linear solver defined! Possible solver block IDs range "
          "from 1-9.");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNewton::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ =
      linesearchfactory.Create(Configuration(), MyGetParameter<std::string>("line search"));

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNewton::ApplyInverse(
    const Epetra_MultiVector& f, Epetra_MultiVector& x) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorNewton::ApplyInverse");
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
  // initialize stuff for Newton loop
  // ---------------------------------------------------------------------------
  // solution increment vector
  Teuchos::RCP<Epetra_MultiVector> inc = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  // residual vector already including negative sign
  Teuchos::RCP<Epetra_MultiVector> fnew = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->ComputeF(x, *fnew);

  // some initializations
  int iter = 0;                                                    // iteration counter
  double steplength = 1.0;                                         // line search parameter
  double fnorm2 = 1.0e+12;                                         // residual L2 norm
  bool converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);  // convergence flag
  bool suffdecr = false;  // flag for sufficient decrease of line search
  bool refactor = false;  // Do we need to re-factor the Jacobian?

  // check for stagnation of iterations
  Teuchos::RCP<NLNSOL::UTILS::StagnationDetection> stagdetect =
      Teuchos::rcp(new NLNSOL::UTILS::StagnationDetection());
  stagdetect->Init(Configuration(),
      MyGetParameter<std::string>("nonlinear operator: stagnation detection"), fnorm2);
  Teuchos::RCP<Teuchos::ParameterList> statusparams = Teuchos::rcp(new Teuchos::ParameterList());

  // print initial state
  PrintIterSummary(iter, fnorm2);

  // ---------------------------------------------------------------------------
  // do a Newton-type iteration loop
  // ---------------------------------------------------------------------------
  while (ContinueIterations(iter, converged))
  {
    ++iter;

    UpdateJacobian(refactor);

    ComputeSearchDirection(fnew, inc, refactor);
    ComputeStepLength(x, *fnew, *inc, fnorm2, steplength, suffdecr);

    // Iterative update
    err = x.Update(steplength, *inc, 1.0);
    if (err != 0)
    {
      dserror("Failed.");
    }

    // evaluate and check for convergence
    NlnProblem()->ComputeF(x, *fnew);
    converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);

    // check for stagnation
    stagdetect->Check(fnorm2);

    stagdetect->Status(statusparams);

    //    if (statusparams->get<double>("Stagnation Detection: ratio") >= 0.9999999999)
    //    {
    //      // revert previous iteraton step
    //      err = x.Update(-steplength, *inc, 1.0);
    //      if (err != 0) { dserror("Failed."); }
    //
    //      // evaluate and check for convergence
    //      NlnProblem()->ComputeF(x, *fnew);
    //      converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);
    //    }

    if (converged)  // or stagdetect->Status() or oparams->get<double>("Stagnation Detection:
                    // ratio") >= 0.9999999999)
    {
      PrintIterSummary(iter, fnorm2);
      break;
    }

    PrintIterSummary(iter, fnorm2);
  }

  bool stagnation = false;
  stagdetect->Status(statusparams);

  if (statusparams->get<bool>("Stagnation Detection: status") == true)
  //      or statusparams->get<double>("Stagnation Detection: ratio") >= 0.9999999999)
  {
    stagnation = true;
  }

  // ---------------------------------------------------------------------------
  // Finish ApplyInverse()
  // ---------------------------------------------------------------------------
  // determine error code
  NLNSOL::UTILS::OperatorStatus errorcode = ErrorCode(iter, converged, err, stagnation);

  if (stagnation) errorcode = NLNSOL::UTILS::opstatus_stagnation;

  // write to output parameter list
  SetOutParameterIter(iter);
  SetOutParameterResidualNorm(fnorm2);
  SetOutParameterConverged(converged);
  SetOutParameterErrorCode(errorcode);
  SetOutParameterStagnation(statusparams);

  // return error code
  return errorcode;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNewton::ComputeSearchDirection(Teuchos::RCP<Epetra_MultiVector> rhs,
    Teuchos::RCP<Epetra_MultiVector> inc, const bool refactor) const
{
  // compute search direction with either fixed or most recent, updated jacobian
  int linsolve_error = linsolver_->Solve(
      NlnProblem()->GetJacobianOperator(), inc, rhs, refactor, refactor, Teuchos::null);

  // test for failure of linear solver
  if (IsSolver() and linsolve_error != 0)
  {
    dserror("Linear solver failed.");
  }

  return linsolve_error;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNewton::ComputeStepLength(const Epetra_MultiVector& x,
    const Epetra_MultiVector& f, const Epetra_MultiVector& inc, double fnorm2, double& lsparam,
    bool& suffdecr) const
{
  const std::string lslist = MyGetParameter<std::string>("line search");

  linesearch_->Init(NlnProblem(), Configuration(), lslist, x, f, inc, fnorm2);
  linesearch_->Setup();
  linesearch_->ComputeLSParam(lsparam, suffdecr);

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNewton::UpdateJacobian(bool& refactor) const
{
  switch (jacevery_)
  {
    case 0:
    {
      refactor = false;
      break;
    }
    case 1:
    {
      NlnProblem()->ComputeJacobian();
      refactor = true;
      break;
    }
    default:
    {
      dserror(
          "Updating the Jacobian every %i iterations is not allowed. "
          "Choose '0' for Chord's method and '1' for Full Newton.",
          jacevery_);
      break;
    }
  }

  return;
}
