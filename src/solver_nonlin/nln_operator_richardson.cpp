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
: linprec_(Teuchos::null)
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

  // ---------------------------------------------------------------------------
  // Get parameter list for linear preconditioner from input file
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> solverparams =
      CreateIFPACKParameterList();


  // ---------------------------------------------------------------------------
  // Create the linear preconditioner based on user's choice
  // ---------------------------------------------------------------------------
  if (solverparams->isSublist("IFPACK Parameters"))
  {
    linprec_ = Teuchos::rcp(
        new LINALG::SOLVER::IFPACKPreconditioner(
            DRT::Problem::Instance()->ErrorFile()->Handle(),
            solverparams->sublist("IFPACK Parameters"),
            solverparams->sublist("Aztec Parameters")));
  }
  else
  {
    dserror("This type of linear preconditioner is not implemented, yet. "
        "Or maybe you just chose a wrong configuration in the input file.");
  }

  Teuchos::RCP<Epetra_MultiVector> tmp =
      Teuchos::rcp(new Epetra_MultiVector(NlnProblem()->DofRowMap(), true));

  // Setup with dummy vectors
  linprec_->Setup(true, &*(NlnProblem()->GetJacobianOperator()),
      &*tmp, &*tmp);

  // Setup() has been called
  SetIsSetup();

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
  double omega =
      Params().get<double>("Richardson iteration: relaxation parameter"); // relaxation parameter

  while (iter < GetMaxIter())
  {
    ++iter;

    err = NlnProblem()->GetJacobianOperator()->Apply(*dx, *r);
    if (err != 0) { dserror("Apply failed."); }

    err = r->Update(-1.0, *fnew, -1.0);
    if (err != 0) { dserror("Update failed."); }

    // compute new iterative increment
    err = linprec_->ApplyInverse(*r, *inc);
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
  dserror("There is no line search algorithm available for a Richardson "
      "iteration preconditioner object.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList>
NLNSOL::NlnOperatorRichardson::CreateIFPACKParameterList() const
{
  // ---------------------------------------------------------------------------
  // Get parameter list for linear preconditioner from input file
  // ---------------------------------------------------------------------------
  // define parameter list (to be filled)
  Teuchos::RCP<Teuchos::ParameterList> params =
      Teuchos::rcp(new Teuchos::ParameterList());

  // get the solver number to read from input file
  const int linsolvernumber = Params().get<int>(
      "Richardson iteration: linear solver");

  // check if the solver ID is valid
  if (linsolvernumber == (-1)) { dserror("No valid linear solver defined!"); }

  // linear solver parameter list
  const Teuchos::ParameterList& inputparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);

  // ---------------------------------------------------------------------------
  // Perform some safety checks on the parameter configuration
  // ---------------------------------------------------------------------------
  const std::string prectype =
      inputparams.get<std::string>("AZPREC");
  if (prectype == "Jacobi" or prectype == "SymmGaussSeidel")
  {
    if (inputparams.get<int>("IFPACKGFILL") != 1)
    {
      dserror("Iterations are performed by %s. "
          "Set input parameter IFPACKGFILL to 1.", Label());
    }
  }

  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Transform Baci solver parameters to Trilinos solver parameters
  // ---------------------------------------------------------------------------
  *params = LINALG::Solver::TranslateSolverParameters(inputparams);

  // ---------------------------------------------------------------------------

  return params;
}
