/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_fas.cpp

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
#include <Epetra_Vector.h>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "fas_hierarchy.H"
#include "fas_nlnlevel.H"
#include "nln_operator_fas.H"
#include "nln_problem.H"
#include "nln_problem_coarselevel.H"
#include "nln_utils.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorFas::NlnOperatorFas()
: hierarchy_(Teuchos::null),
  cycletype_(INPAR::NLNSOL::FAS::cycle_none)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorFas::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnOperatorFas::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // create the multigrid level hierarchy
  hierarchy_ = Teuchos::rcp(new NLNSOL::FAS::AMGHierarchy());
  hierarchy_->Init(Comm(), Params(), NlnProblem(),
      BaciLinearSolver()->Params().sublist("ML Parameters"));
  hierarchy_->Setup();

  std::string cycletype = "V"; // ToDo (mayr) read from parameter list
//      Params().sublist("FAS: MueLu Parameters").get<std::string>("cycle type");
  if (cycletype == "V")
    cycletype_ = INPAR::NLNSOL::FAS::cycle_v;
  else if (cycletype == "W")
    cycletype_ = INPAR::NLNSOL::FAS::cycle_w;
  else
    dserror("Unknown multigrid cycle type '%s'", cycletype.c_str());

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorFas::RebuildPrecConst() const
{
  if (hierarchy_.is_null())
    dserror("Hierarchy of multigrid levels 'hierarchy_' is not set, yet.");

  NlnProblem()->ComputeJacobian();
  hierarchy_->RefreshRAPs();

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::ApplyInverse(const Epetra_MultiVector& f_do_not_use,
    Epetra_MultiVector& x) const
{
  int err = 0;

  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnOperatorFas::ApplyInverse");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() and Setup() have been called
  if (not IsInit())  { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // local copy of residual that can be modified
  Teuchos::RCP<Epetra_MultiVector> f =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->ComputeF(x, *f);

  // local copy of solution to work with
  Teuchos::RCP<Epetra_MultiVector> xtmp =
      Teuchos::rcp(new Epetra_MultiVector(x));

  double fnorm2 = 1.0e+12;
  bool converged = NlnProblem()->ConvergenceCheck(*f, fnorm2);

  *getOStream() << std::setprecision(12) << std::scientific;
  *getOStream() << "Starting FAS with residual norm = " << fnorm2 << std::endl;

  Teuchos::RCP<NLNSOL::UTILS::StagnationDetection> stagdetect =
      Teuchos::rcp(new NLNSOL::UTILS::StagnationDetection());
  stagdetect->Init(
      Params().sublist("Nonlinear Operator: Stagnation Detection"), fnorm2);

  int iter = 0;
  while (ContinueIterations(iter, converged))
  {
    ++iter;

    if (getVerbLevel() > Teuchos::VERB_NONE)
    {
      *getOStream() << "Start multigrid cycle for the " << iter << ". time."
          << std::endl;
    }

    // call generic cycling routine
    Cycle(xtmp, 0);

    // Evaluate and check for convergence
    NlnProblem()->ComputeF(*xtmp, *f);
    converged = NlnProblem()->ConvergenceCheck(*f, fnorm2);

    stagdetect->Check(fnorm2);

    RebuildPrecConst();

    PrintIterSummary(iter, fnorm2);
  }

  // Update reference to solution
  err = x.Update(1.0, *xtmp, 0.0);
  if (err != 0) { dserror("Update failed."); }

  bool stagnation = false;
  if (stagdetect->Status() or Hierarchy()->CheckAllLevelStagnation())
    stagnation = true;

  // ---------------------------------------------------------------------------
  // Finish ApplyInverse()
  // ---------------------------------------------------------------------------
  // determine error code
  NLNSOL::UTILS::OperatorStatus errorcode =
      ErrorCode(iter, converged, stagnation);

  // write to output parameter list
  SetOutParameterIter(iter);
  SetOutParameterResidualNorm(fnorm2);
  SetOutParameterConverged(converged);
  SetOutParameterErrorCode(errorcode);

  Teuchos::RCP<Teuchos::ParameterList> status =
      Teuchos::rcp(new Teuchos::ParameterList(*stagdetect->StatusParams()));
  status->set<bool>("Stagnation Detection: status", stagnation);
  SetOutParameterStagnation(status);

  if (converged)
    *getOStream() << Label() << " seems to be converged in " << iter
        << " iterations." << std::endl;

  // return error code
  return errorcode;
}

/*----------------------------------------------------------------------------*/
const int NLNSOL::NlnOperatorFas::Cycle(Teuchos::RCP<Epetra_MultiVector>& x,
    const int level) const
{
  // error code
  int err = 0;

  // chose multigrid cycle type based on user input
  switch (cycletype_)
  {
  case INPAR::NLNSOL::FAS::cycle_none:
  {
    dserror("No multigrid cycle type chosen. Choose a valid cycle type!");
    err = -1;
    break;
  }
  case INPAR::NLNSOL::FAS::cycle_v:
  {
    err = VCycle(x, level);
    break;
  }
  case INPAR::NLNSOL::FAS::cycle_w:
  {
    err = WCycle();
    break;
  }
  default:
  {
    dserror("Unknown multigrid cycle type. Choose a valid cycle type!");
    err = -1;
    break;
  }
  }

  return err;
}

/*----------------------------------------------------------------------------*/
const int NLNSOL::NlnOperatorFas::VCycle(Teuchos::RCP<Epetra_MultiVector>& xbar,
    const int level) const
{
  *getOStream() << "Entering VCycle on level " << level << std::endl;

  // error code
  int err = 0;

  Hierarchy()->CheckLevelID(level); // check for valid level ID

  /* Copy solution to work with on this level, since 'xbar' has to be kept
   * untouched in order to use it later to compute the coarse grid correction.
   */
  Teuchos::RCP<Epetra_MultiVector> x =
      Teuchos::rcp(new Epetra_MultiVector(*xbar));

  if (level < Hierarchy()->NumLevels() - 1)
  {
    PreSmoothing(*x, level);

    // recursive call of V-cycle
    {
      // -----------------------------------------------------------------------
      // prepare next coarser level
      // -----------------------------------------------------------------------
      // restrict current solution (this will be 'xbar' on next coarser level)
      Teuchos::RCP<Epetra_MultiVector> xc =
          Teuchos::rcp(new Epetra_MultiVector(*x));
      xc = Hierarchy()->NlnLevel(level + 1)->RestrictToNextCoarserLevel(xc);

      // restrict current residual (this will be 'fbar' on next coarser level)
      Teuchos::RCP<Epetra_MultiVector> fc =
          Hierarchy()->NlnLevel(level)->NlnProblem()->ComputePlainF(*x);
      fc = Hierarchy()->NlnLevel(level + 1)->RestrictToNextCoarserLevel(fc);

      Hierarchy()->NlnLevel(level + 1)->NlnProblem()->SetupResidualModification(
          xc, fc);
      // -----------------------------------------------------------------------

      err = VCycle(xc, level + 1);
      if (err != 0)
        dserror("VCycle on level %d returned error %d.", level+1, err);

      // apply coarse level correction
      xc = Hierarchy()->NlnLevel(level+1)->ProlongateToNextFinerLevel(xc);
      err = x->Update(1.0, *xc, 1.0);
      if (err != 0) { dserror("Update failed."); }
    }

    PostSmoothing(*x, level);
  }
  else
  {
    CoarseLevelSolve(*x, level);
  }

  // update solution variable  'xbar' (to be returned)
  if (level > 0) // compute coarse level correction
  {
    err = xbar->Update(1.0, *x, -1.0);
    if (err != 0) { dserror("Update failed."); }
  }
  else // fine level: copy solution to return variable 'xbar'
  {
    err = xbar->Update(1.0, *x, 0.0);
    if (err != 0) { dserror("Update failed."); }
  }

  *getOStream() << "Leaving VCycle on level " << level << std::endl;

  return err;
}

/*----------------------------------------------------------------------------*/
const int NLNSOL::NlnOperatorFas::WCycle() const
{
  dserror("W-Cycle not implemented, yet.");

  return -1;
}

/*----------------------------------------------------------------------------*/
const int NLNSOL::NlnOperatorFas::PreSmoothing(
    Epetra_MultiVector& x, const int level) const
{
  // evaluate current residual
  Teuchos::RCP<Epetra_MultiVector> fsmoothed =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(x, *fsmoothed);

  // pre-smoothing
  return Hierarchy()->NlnLevel(level)->DoPreSmoothing(*fsmoothed, x);
}

/*----------------------------------------------------------------------------*/
const int NLNSOL::NlnOperatorFas::PostSmoothing(
    Epetra_MultiVector& x, const int level) const
{
  // evaluate current residual
  Teuchos::RCP<Epetra_MultiVector> fsmoothed =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(x, *fsmoothed);

  // post-smoothing
  return Hierarchy()->NlnLevel(level)->DoPostSmoothing(*fsmoothed, x);
}

/*----------------------------------------------------------------------------*/
const int NLNSOL::NlnOperatorFas::CoarseLevelSolve(
    Epetra_MultiVector& x, const int level) const
{
  if (getVerbLevel() > Teuchos::VERB_NONE)
  {
      *getOStream() << std::endl
          << "**************************" << std::endl
          << "*** COARSE LEVEL SOLVE ***" << std::endl
          << "**************************" << std::endl
          << std::endl;
  }

  // evaluate current residual // ToDo Do we really need to ComputeF() here?
  Teuchos::RCP<Epetra_MultiVector> f =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(x, *f);
  return Hierarchy()->NlnLevel(level)->DoCoarseLevelSolve(*f, x);
}

/*----------------------------------------------------------------------------*/
const Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy>
NLNSOL::NlnOperatorFas::Hierarchy() const
{
  if (hierarchy_.is_null())
    dserror("Hierarchy of multigrid levels 'hierarchy_' is not set, yet.");

  return hierarchy_;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorFas::RebuildPrec()
{
  if (hierarchy_.is_null())
    dserror("Hierarchy of multigrid levels 'hierarchy_' is not set, yet.");

  NlnProblem()->ComputeJacobian();
  hierarchy_->RefreshRAPs();

  return;
}
