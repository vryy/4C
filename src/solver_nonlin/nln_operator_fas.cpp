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
    : hierarchy_(Teuchos::null), cycletype_(INPAR::NLNSOL::FAS::cycle_none)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorFas::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorFas::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  // create the multigrid level hierarchy
  hierarchy_ = Teuchos::rcp(new NLNSOL::FAS::AMGHierarchy());
  hierarchy_->Init(Comm(), Configuration(), MyGetParameter<std::string>("FAS: AMG Hierarchy"),
      NlnProblem(), BaciLinearSolver()->Params().sublist("ML Parameters"));
  hierarchy_->Setup();

  const std::string cycletype = MyGetParameter<std::string>("FAS: cycle type");
  if (cycletype == "V")
    cycletype_ = INPAR::NLNSOL::FAS::cycle_v;
  else if (cycletype == "V 2-level")
    cycletype_ = INPAR::NLNSOL::FAS::cycle_v_two_level;
  else if (cycletype == "V 3-level")
    cycletype_ = INPAR::NLNSOL::FAS::cycle_v_three_level;
  else if (cycletype == "W")
    cycletype_ = INPAR::NLNSOL::FAS::cycle_w;
  else
    dserror("Unknown multigrid cycle type '%s'", cycletype.c_str());

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::ApplyInverse(
    const Epetra_MultiVector& f_do_not_use, Epetra_MultiVector& x) const
{
  int err = 0;

  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorFas::ApplyInverse");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() and Setup() have been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  // local copy of solution to work with
  Teuchos::RCP<Epetra_MultiVector> xtmp = Teuchos::rcp(new Epetra_MultiVector(x));

  // residual
  Teuchos::RCP<Epetra_MultiVector> f = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->ComputeF(*xtmp, *f);

  // initial residual norm
  double fnorm2 = 1.0e+12;
  bool converged = NlnProblem()->ConvergenceCheck(*f, fnorm2);
  PrintIterSummary(0, fnorm2);

  Teuchos::RCP<NLNSOL::UTILS::StagnationDetection> stagdetect =
      Teuchos::rcp(new NLNSOL::UTILS::StagnationDetection());
  stagdetect->Init(Configuration(),
      MyGetParameter<std::string>("nonlinear operator: stagnation detection"), fnorm2);

  int iter = 0;
  while (ContinueIterations(iter, converged))
  {
    ++iter;

    if (getVerbLevel() > Teuchos::VERB_LOW)
    {
      *getOStream() << "Start multigrid cycle for the " << iter << ". time." << std::endl;
    }

    // call generic cycling routine
    Cycle(xtmp);

    // Evaluate and check for convergence
    NlnProblem()->ComputeF(*xtmp, *f);
    converged = NlnProblem()->ConvergenceCheck(*f, fnorm2);

    stagdetect->Check(fnorm2);

    //    RebuildPrecConst();

    PrintIterSummary(iter, fnorm2);
  }

  // Update reference to solution
  err = x.Update(1.0, *xtmp, 0.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  bool stagnation = false;
  //  if (stagdetect->Status() or Hierarchy()->CheckAllLevelStagnation())
  //    stagnation = true;

  // ---------------------------------------------------------------------------
  // Finish ApplyInverse()
  // ---------------------------------------------------------------------------
  // determine error code
  NLNSOL::UTILS::OperatorStatus errorcode = ErrorCode(iter, converged, stagnation);

  // write to output parameter list
  SetOutParameterIter(iter);
  SetOutParameterResidualNorm(fnorm2);
  SetOutParameterConverged(converged);
  SetOutParameterErrorCode(errorcode);

  Teuchos::RCP<Teuchos::ParameterList> status =
      Teuchos::rcp(new Teuchos::ParameterList(*stagdetect->StatusParams()));
  status->set<bool>("Stagnation Detection: status", stagnation);
  SetOutParameterStagnation(status);

  if (converged and getVerbLevel() > Teuchos::VERB_LOW)
    *getOStream() << Label() << " seems to be converged in " << iter << " iterations." << std::endl;

  // return error code
  return errorcode;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::Cycle(Teuchos::RCP<Epetra_MultiVector> x) const
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
      err = VCycle(x, 0);  // start recursive V-Cycle on level 0
      break;
    }
    case INPAR::NLNSOL::FAS::cycle_v_two_level:
    {
      err = VCycleTwoLevel(x);
      break;
    }
    case INPAR::NLNSOL::FAS::cycle_v_three_level:
    {
      err = VCycleThreeLevel(x);
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
int NLNSOL::NlnOperatorFas::VCycle(Teuchos::RCP<Epetra_MultiVector> xbar, const int level) const
{
  if (getVerbLevel() > Teuchos::VERB_LOW)
    *getOStream() << "Entering recursive V-Cycle on level " << level << std::endl;

  // error code
  int err = 0;

  Hierarchy()->CheckLevelID(level);  // check for valid level ID

  /* Copy solution to work with on this level, since 'xbar' has to be kept
   * untouched in order to use it later to compute the coarse grid correction.
   */
  Teuchos::RCP<Epetra_MultiVector> x = Teuchos::rcp(new Epetra_MultiVector(*xbar));

  if (level < Hierarchy()->NumLevels() - 1)
  {
    PreSmoothing(*x, level);

    // recursive call of V-cycle
    {
      // -----------------------------------------------------------------------
      // prepare next coarser level
      // -----------------------------------------------------------------------
      // restrict current solution (this will be 'xbar' on next coarser level)
      Teuchos::RCP<Epetra_MultiVector> xc = Teuchos::rcp(new Epetra_MultiVector(*x));
      xc = Hierarchy()->NlnLevel(level + 1)->RestrictToNextCoarserLevel(xc);

      // restrict current residual (this will be 'fbar' on next coarser level)
      Teuchos::RCP<const Epetra_MultiVector> fc =
          Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(*x);
      fc = Hierarchy()->NlnLevel(level + 1)->RestrictToNextCoarserLevel(fc);

      Hierarchy()->NlnLevel(level + 1)->NlnProblem()->SetupResidualModification(xc, fc);
      // -----------------------------------------------------------------------

      err = VCycle(xc, level + 1);
      if (err != 0) dserror("VCycle on level %d returned error %d.", level + 1, err);

      // apply coarse level correction
      xc = Hierarchy()->NlnLevel(level + 1)->ProlongateToNextFinerLevel(xc);
      err = x->Update(1.0, *xc, 1.0);
      if (err != 0)
      {
        dserror("Update failed.");
      }
    }

    PostSmoothing(*x, level);
  }
  else
  {
    CoarseLevelSolve(*x, level);
  }

  // update solution variable  'xbar' (to be returned)
  if (level > 0)  // compute coarse level correction
  {
    err = xbar->Update(1.0, *x, -1.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }
  }
  else  // fine level: copy solution to return variable 'xbar'
  {
    err = xbar->Update(1.0, *x, 0.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }
  }

  if (getVerbLevel() > Teuchos::VERB_LOW)
    *getOStream() << "Leaving VCycle on level " << level << std::endl;

  return err;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::VCycleTwoLevel(Teuchos::RCP<Epetra_MultiVector> x) const
{
  if (getVerbLevel() > Teuchos::VERB_LOW)
    *getOStream() << "Starting explicit 2-level V-Cycle." << std::endl;

  // error code
  int err = 0;

  // level 0
  {
    PreSmoothing(*x, 0);
  }

  // level 1 (= coarse level)
  Teuchos::RCP<Epetra_MultiVector> xc =  // solution vector on coarse level
      Teuchos::rcp(new Epetra_MultiVector(*x));
  {
    // -----------------------------------------------------------------------
    // prepare next coarser level
    // -----------------------------------------------------------------------
    // restrict current solution
    xc = Hierarchy()->NlnLevel(1)->RestrictToNextCoarserLevel(xc);

    // store 'xc' as 'xbar'
    Teuchos::RCP<const Epetra_MultiVector> xbar = Teuchos::rcp(new Epetra_MultiVector(*xc));

    // restrict current residual (this will be 'fbar' on next coarser level)
    Teuchos::RCP<const Epetra_MultiVector> fc =
        Hierarchy()->NlnLevel(0)->NlnProblem()->ComputeF(*x);
    fc = Hierarchy()->NlnLevel(1)->RestrictToNextCoarserLevel(fc);

    Hierarchy()->NlnLevel(1)->NlnProblem()->SetupResidualModification(xc, fc);

    CoarseLevelSolve(*xc, 1);

    // compute coarse level correction
    err = xc->Update(-1.0, *xbar, 1.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }
  }

  // level 0
  {
    // apply coarse level correction
    xc = Hierarchy()->NlnLevel(1)->ProlongateToNextFinerLevel(xc);
    err = x->Update(1.0, *xc, 1.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }

    PostSmoothing(*x, 0);
  }

  if (getVerbLevel() > Teuchos::VERB_LOW)
    *getOStream() << "Leaving explicit 2-level V-Cycle." << std::endl;

  return 0;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::VCycleThreeLevel(Teuchos::RCP<Epetra_MultiVector> x) const
{
  if (getVerbLevel() > Teuchos::VERB_LOW)
    *getOStream() << "Starting explicit 3-level V-Cycle." << std::endl;

  // error code
  int err = 0;

  // level 0 (= fine level)
  PreSmoothing(*x, 0);

  // level 1 (= medium level)
  Teuchos::RCP<Epetra_MultiVector> xc =  // solution vector on medium level
      Teuchos::rcp(new Epetra_MultiVector(*x));
  // -----------------------------------------------------------------------
  // prepare next coarser level
  // -----------------------------------------------------------------------
  // restrict current solution
  xc = Hierarchy()->NlnLevel(1)->RestrictToNextCoarserLevel(xc);

  // store 'xc' as 'xbar'
  Teuchos::RCP<const Epetra_MultiVector> xcbar = Teuchos::rcp(new Epetra_MultiVector(*xc));

  // restrict current residual (this will be 'fbar' on next coarser level)
  Teuchos::RCP<const Epetra_MultiVector> fc = Hierarchy()->NlnLevel(0)->NlnProblem()->ComputeF(*x);
  fc = Hierarchy()->NlnLevel(1)->RestrictToNextCoarserLevel(fc);

  Hierarchy()->NlnLevel(1)->NlnProblem()->SetupResidualModification(xc, fc);

  PreSmoothing(*xc, 1);

  // level 2 (= coarse level)
  Teuchos::RCP<Epetra_MultiVector> xcc =  // solution vector on coarse level
      Teuchos::rcp(new Epetra_MultiVector(*xc));
  // -----------------------------------------------------------------------
  // prepare next coarser level
  // -----------------------------------------------------------------------
  // restrict current solution
  xcc = Hierarchy()->NlnLevel(2)->RestrictToNextCoarserLevel(xcc);

  // store 'xc' as 'xbar'
  Teuchos::RCP<const Epetra_MultiVector> xccbar = Teuchos::rcp(new Epetra_MultiVector(*xcc));

  // restrict current residual (this will be 'fbar' on next coarser level)
  Teuchos::RCP<const Epetra_MultiVector> fcc =
      Hierarchy()->NlnLevel(1)->NlnProblem()->ComputeF(*xc);
  fcc = Hierarchy()->NlnLevel(2)->RestrictToNextCoarserLevel(fcc);

  Hierarchy()->NlnLevel(2)->NlnProblem()->SetupResidualModification(xcc, fcc);

  CoarseLevelSolve(*xcc, 2);

  // compute coarse level correction
  err = xcc->Update(-1.0, *xccbar, 1.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  // apply coarse level correction
  xcc = Hierarchy()->NlnLevel(2)->ProlongateToNextFinerLevel(xcc);
  err = xc->Update(1.0, *xcc, 1.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  PostSmoothing(*xc, 1);

  // compute coarse level correction
  err = xc->Update(-1.0, *xcbar, 1.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  // level 0 (= fine level)
  // apply coarse level correction
  xc = Hierarchy()->NlnLevel(1)->ProlongateToNextFinerLevel(xc);
  err = x->Update(1.0, *xc, 1.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  PostSmoothing(*x, 0);

  if (getVerbLevel() > Teuchos::VERB_LOW)
    *getOStream() << "Leaving explicit 3-level V-Cycle." << std::endl;

  return 0;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::WCycle() const
{
  dserror("W-Cycle not implemented, yet.");

  return -1;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::PreSmoothing(Epetra_MultiVector& x, const int level) const
{
  // evaluate current residual
  Teuchos::RCP<Epetra_MultiVector> fsmoothed = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(x, *fsmoothed);

  // pre-smoothing
  return Hierarchy()->NlnLevel(level)->DoPreSmoothing(*fsmoothed, x);
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::PostSmoothing(Epetra_MultiVector& x, const int level) const
{
  // evaluate current residual
  Teuchos::RCP<Epetra_MultiVector> fsmoothed = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(x, *fsmoothed);

  // post-smoothing
  return Hierarchy()->NlnLevel(level)->DoPostSmoothing(*fsmoothed, x);
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorFas::CoarseLevelSolve(Epetra_MultiVector& x, const int level) const
{
  if (getVerbLevel() > Teuchos::VERB_LOW)
  {
    *getOStream() << std::endl
                  << "**************************" << std::endl
                  << "*** COARSE LEVEL SOLVE ***" << std::endl
                  << "**************************" << std::endl
                  << std::endl;
  }

  // evaluate current residual // ToDo Do we really need to ComputeF() here?
  Teuchos::RCP<Epetra_MultiVector> f = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  Hierarchy()->NlnLevel(level)->NlnProblem()->ComputeF(x, *f);
  return Hierarchy()->NlnLevel(level)->DoCoarseLevelSolve(*f, x);
}

/*----------------------------------------------------------------------------*/
const Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> NLNSOL::NlnOperatorFas::Hierarchy() const
{
  if (hierarchy_.is_null()) dserror("Hierarchy of multigrid levels 'hierarchy_' is not set, yet.");

  return hierarchy_;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorFas::RebuildPrec()
{
  if (hierarchy_.is_null()) dserror("Hierarchy of multigrid levels 'hierarchy_' is not set, yet.");

  NlnProblem()->ComputeJacobian();
  hierarchy_->RefreshRAPs();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorFas::RebuildPrecConst() const
{
  if (hierarchy_.is_null()) dserror("Hierarchy of multigrid levels 'hierarchy_' is not set, yet.");

  NlnProblem()->ComputeJacobian();
  hierarchy_->RefreshRAPs();

  return;
}
