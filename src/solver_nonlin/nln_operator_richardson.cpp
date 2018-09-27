/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_richardson.cpp

\brief Nonlinear preconditioned Richardson iterations with damping

\level 3

\maintainer Matthias Mayr
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

// Ifpack
#include <Ifpack.h>

// standard
#include <iostream>
#include <vector>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "nln_operator_factory.H"
#include "nln_operator_richardson.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"

#include "../solver/solver_preconditionertype.H"
#include "../solver/solver_ifpackpreconditioner.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorRichardson::NlnOperatorRichardson() : nlnprec_(Teuchos::null) { return; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorRichardson::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorRichardson::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  SetupPreconditioner();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorRichardson::SetupPreconditioner()
{
  const std::string opname =
      MyGetParameter<std::string>("nonlinear richardson: nonlinear preconditioner");

  NlnOperatorFactory nlnopfactory;
  nlnprec_ = nlnopfactory.Create(Configuration(), opname);
  nlnprec_->Init(Comm(), Configuration(), opname, NlnProblem(), BaciLinearSolver(), Nested() + 1);
  nlnprec_->Setup();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorRichardson::ApplyInverse(
    const Epetra_MultiVector& f, Epetra_MultiVector& x) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorRichardson::ApplyInverse");
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

  Teuchos::RCP<Epetra_MultiVector> xprec = Teuchos::rcp(new Epetra_MultiVector(x));
  Teuchos::RCP<Epetra_MultiVector> inc = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  // initial residual
  Teuchos::RCP<Epetra_MultiVector> fnew = Teuchos::rcp(new Epetra_MultiVector(xprec->Map(), true));
  NlnProblem()->ComputeF(*xprec, *fnew);
  double fnorm2 = 1.0e+12;
  bool converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);
  PrintIterSummary(0, fnorm2);

  Teuchos::RCP<NLNSOL::UTILS::StagnationDetection> stagdetect =
      Teuchos::rcp(new NLNSOL::UTILS::StagnationDetection());
  stagdetect->Init(Configuration(),
      MyGetParameter<std::string>("nonlinear operator: stagnation detection"), fnorm2);

  // damping parameter
  const double omega = MyGetParameter<double>("nonlinear richardson: damping parameter");

  int iter = 0;
  while (ContinueIterations(iter, converged))
  {
    ++iter;

    nlnprec_->ApplyInverse(*fnew, *xprec);
    err = inc->Update(1.0, *xprec, -1.0, x, 0.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }

    err = x.Update(omega, *inc, 1.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }
    err = xprec->Update(1.0, x, 0.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }

    NlnProblem()->ComputeF(x, *fnew);
    converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);
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
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorRichardson::ComputeStepLength(const Epetra_MultiVector& x,
    const Epetra_MultiVector& f, const Epetra_MultiVector& inc, double fnorm2, double& lsparam,
    bool& suffdecr) const
{
  if (MyGetParameter<std::string>("nonlinear richardson: damping strategy") ==
      "fixed damping parameter")
  {
    lsparam = MyGetParameter<double>("nonlinear richardson: damping parameter");
  }
  else if (MyGetParameter<std::string>("nonlinear richardson: damping strategy") == "line search")
    dserror("Not implemented, yet.");
  else
    dserror("Unknown strategy ");

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorRichardson::RebuildPrec()
{
  nlnprec_->RebuildPrec();

  return;
}
