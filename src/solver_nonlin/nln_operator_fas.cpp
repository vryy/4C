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
  hierarchy_->Init(Comm(), Params(), NlnProblem());
  hierarchy_->Setup();

  std::string cycletype = Params().sublist("FAS: MueLu Parameters").get<std::string>("cycle type");
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
int NLNSOL::NlnOperatorFas::ApplyInverse(const Epetra_MultiVector& f,
    Epetra_MultiVector& x) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnOperatorFas::ApplyInverse");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() and Setup() have been called
  if (not IsInit())  { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // local copy of residual that can be modified
  Teuchos::RCP<Epetra_MultiVector> ftmp =
      Teuchos::rcp(new Epetra_MultiVector(f));

  bool converged = false;
  double fnorm2 = 1.0e+12;
  NlnProblem()->ConvergenceCheck(f, fnorm2);

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
    Cycle(*ftmp, x, 0);

    // Evaluate and check for convergence
    NlnProblem()->ComputeF(x, *ftmp);
    converged = NlnProblem()->ConvergenceCheck(*ftmp, fnorm2);

    stagdetect->Check(fnorm2);

    PrintIterSummary(iter, fnorm2);
  }

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

  // return error code
  return errorcode;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorFas::Cycle(const Epetra_MultiVector& f,
    Epetra_MultiVector& x,
    const int level
    ) const
{
  // chose multigrid cycle type based on user input
  switch (cycletype_)
  {
  case INPAR::NLNSOL::FAS::cycle_none:
  {
    dserror("No multigrid cycle type chose.");
    break;
  }
  case INPAR::NLNSOL::FAS::cycle_v:
  {
    VCycle(f, x, level);
    break;
  }
  case INPAR::NLNSOL::FAS::cycle_w:
  {
    WCycle();
    break;
  }
  default:
  {
    dserror("Unknown multigrid cycle type.");
    break;
  }
  }

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorFas::VCycle(const Epetra_MultiVector& f,
    Epetra_MultiVector& x,
    const int level
    ) const
{
  int err = 0;

  Hierarchy()->CheckLevelID(level);

  if (getVerbLevel() > Teuchos::VERB_NONE)
    *getOStream() << "WELCOME to VCycle on level " << level << "." << std::endl;

  // we need at least zeroed vectors, especially on the fine level
  // restriction of fine-level residual
  Teuchos::RCP<Epetra_MultiVector> fbar =
      Teuchos::rcp(new Epetra_MultiVector(f));
  // restriction of fine-level solution
  Teuchos::RCP<Epetra_MultiVector> xbar =
      Teuchos::rcp(new Epetra_MultiVector(x));
  // coarse-grid evaluation of residual
  Teuchos::RCP<Epetra_MultiVector> fhat =
      Teuchos::rcp(new Epetra_MultiVector(f.Map(), 1, true));

  /* Leave 'x' untouched on the level since it is needed to approximate the
   * error after the post-smoothing. Use 'xtemp' instead for all operations on
   * this level.
   */
  Teuchos::RCP<Epetra_MultiVector> xtemp =
      Teuchos::rcp(new Epetra_MultiVector(x));

  // compute coarse level contributions
  if (level > 0)
  {
    // restrict fine level residual to current (coarser) level
    Hierarchy()->NlnLevel(level)->RestrictToNextCoarserLevel(fbar);

    // restrict fine level solution to current (coarser) level
    Hierarchy()->NlnLevel(level)->RestrictToNextCoarserLevel(xbar);
    Hierarchy()->NlnLevel(level)->RestrictToNextCoarserLevel(fhat);

    // Set zero vectors for fhat and fbar to enable evaluation of fhat.
    // Note: 'fhat' is zero here.
    Hierarchy()->NlnLevel(level)->NlnProblem()->SetFHatFBar(fhat, fhat);

    // evaluate coarse-grid residual
    Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(*xbar, *fhat);

    Hierarchy()->NlnLevel(level)->RestrictToNextCoarserLevel(xtemp);

    // set fhat_ and fbar_ in current level
    Hierarchy()->NlnLevel(level)->NlnProblem()->SetFHatFBar(fhat, fbar);
  }
  else if (level == 0) // do nothing on the fine level
  {
    fbar->PutScalar(0.0);
    xbar->PutScalar(0.0);
    fhat->PutScalar(0.0);
  } // do nothing on the fine level
  else
  {
    dserror("Level ID %d is not a valid level ID. Level ID has to be in [0,%d]",
        level, Hierarchy()->NumLevels() - 1);
  }

#ifdef DEBUG
  // additional safety checks w.r.t. maps
  if (not Hierarchy()->NlnLevel(level)->DofRowMap().PointSameAs(fbar->Map()))
    dserror("Maps do not match.");
  if (not Hierarchy()->NlnLevel(level)->DofRowMap().PointSameAs(xbar->Map()))
    dserror("Maps do not match.");
  if (not Hierarchy()->NlnLevel(level)->DofRowMap().PointSameAs(fhat->Map()))
    dserror("Maps do not match.");
  if (not Hierarchy()->NlnLevel(level)->DofRowMap().PointSameAs(xtemp->Map()))
    dserror("Maps do not match.");
#endif

  // do further coarsening only in case we are not on the coarsest level, yet.
  if (level + 1 < Hierarchy()->NumLevels())
  {
    // evaluate current residual
    Teuchos::RCP<Epetra_MultiVector> fsmoothed =
        Teuchos::rcp(new Epetra_MultiVector(xtemp->Map(), true));
    Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(*xtemp, *fsmoothed);

    // pre-smoothing
    Hierarchy()->NlnLevel(level)->DoPreSmoothing(*fsmoothed, *xtemp);

    // evaluate current residual
    Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(*xtemp, *fsmoothed);

    // call VCycle on next coarser level recursively
    VCycle(*fsmoothed, *xtemp, level + 1);

    // evaluate current residual
    Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(*xtemp, *fsmoothed);

    // post-smoothing
    Hierarchy()->NlnLevel(level)->DoPostSmoothing(*fsmoothed, *xtemp);
  }
  else // coarse level solve
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
    Teuchos::RCP<Epetra_MultiVector> fsmoothed =
        Teuchos::rcp(new Epetra_MultiVector(xtemp->Map(), true));
    Hierarchy()->NlnLevel(Hierarchy()->NumLevels() - 1)->NlnProblem()->ComputeF(
        *xtemp, *fsmoothed);
    Hierarchy()->NlnLevel(Hierarchy()->NumLevels() - 1)->DoCoarseLevelSolve(
        *fsmoothed, *xtemp);
  }

  if (level > 0)
  {
    // coarse-grid correction (=approximation of error)
    Teuchos::RCP<Epetra_MultiVector> correction =
        Teuchos::rcp(new Epetra_MultiVector(xtemp->Map(), 1, true));
    err = correction->Update(1.0, *xtemp, -1.0, *xbar, 0.0);
    if (err != 0) { dserror("Failed!"); }

    // prolongate correction to finer level and apply it
    err = Hierarchy()->NlnLevel(level)->ProlongateToNextFinerLevel(correction);
    if (err != 0) { dserror("Failed!"); }
    x.Update(1.0, *correction, 1.0);
  }
  else /* no correction on the finest level */
  {
    x.Update(1.0, *xtemp, 0.0);
  }

#ifdef DEBUG
  // check whether all level transfers handled the maps correctly
  if (not x.Map().PointSameAs(f.Map()))
    dserror("Map failure during recursive calls of V-cycle!");
#endif

  if (getVerbLevel() > Teuchos::VERB_NONE)
  {
      *getOStream() << "GOOD BYE from VCycle on level " << level << "."
          <<  std::endl;
  }

  return;
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
