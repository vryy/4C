/*----------------------------------------------------------------------------*/
/*!
 \file nln_operator_richardson.cpp

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
NLNSOL::NlnOperatorRichardson::NlnOperatorRichardson()
: nlnprec_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorRichardson::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnOperatorRichardson::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  SetupPreconditioner();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorRichardson::SetupPreconditioner()
{
  const std::string opname = MyGetParameter<std::string>(
      "nonlinear richardson: nonlinear preconditioner");

  NlnOperatorFactory nlnopfactory;
  nlnprec_ = nlnopfactory.Create(Configuration(), opname);
  nlnprec_->Init(Comm(), Configuration(), opname, NlnProblem(),
      BaciLinearSolver(), Nested() + 1);
  nlnprec_->Setup();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorRichardson::ApplyInverse(const Epetra_MultiVector& f,
    Epetra_MultiVector& x) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnOperatorRichardson::ApplyInverse");
  Teuchos::TimeMonitor monitor(*time);

  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // ---------------------------------------------------------------------------
  // initialize stuff for Richardson iteration loop
  // ---------------------------------------------------------------------------
  // increment to solve for
  Teuchos::RCP<Epetra_MultiVector> dx =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  // evaluate at current solution
  Teuchos::RCP<Epetra_MultiVector> fnew =
      Teuchos::rcp(new Epetra_MultiVector(f.Map(), true));
  NlnProblem()->ComputeF(x, *fnew);

  // linear residual
  Teuchos::RCP<Epetra_MultiVector> r =
      Teuchos::rcp(new Epetra_MultiVector(*fnew));

  // auxiliary iterative increment to feed to preconditioner
  Teuchos::RCP<Epetra_MultiVector> inc =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  int iter = 0; // iteration counter
  double omega = MyGetParameter<double>(
      "nonlinear richardson: relaxation parameter"); // relaxation parameter

  while (iter < GetMaxIter())
  {
    ++iter;

    err = NlnProblem()->GetJacobianOperator()->Apply(*dx, *r);
    if (err != 0) { dserror("Apply failed."); }

    err = r->Update(-1.0, *fnew, -1.0);
    if (err != 0) { dserror("Update failed."); }

    // compute new iterative increment
    err = nlnprec_->ApplyInverse(*r, *inc);
    if (err != 0) { dserror("ApplyInverse() failed."); }

    // update solution increment
    err = dx->Update(omega, *inc, 1.0);
    if (err != 0) { dserror("Update failed."); }
  }

  // update solution
  err = x.Update(1.0, *dx, 1.0);
  if (err != 0) { dserror("Update failed."); }

  // evaluate at the new solution
  NlnProblem()->ComputeF(x, *fnew);
  double fnorm2 = 1.0e+12;
  bool converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);
  PrintIterSummary(iter, fnorm2);

  // ---------------------------------------------------------------------------
  // Finish ApplyInverse()
  // ---------------------------------------------------------------------------
  // determine error code
  NLNSOL::UTILS::OperatorStatus errorcode =
      ErrorCode(iter, converged, err);

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
    const Epetra_MultiVector& f, const Epetra_MultiVector& inc, double fnorm2,
    double& lsparam, bool& suffdecr) const
{
  dserror("There is no line search algorithm available for nonlinear "
      "Richardson iteration.");

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorRichardson::RebuildPrec()
{
  nlnprec_->RebuildPrec();

  return;
}
