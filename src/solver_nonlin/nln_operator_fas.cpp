/*----------------------------------------------------------------------*/
/*!
\file nln_operator_fas.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "fas_hierarchy.H"
#include "fas_nlnlevel.H"
#include "nln_operator_fas.H"
#include "nln_problem.H"
#include "nln_problem_coarselevel.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperatorFas::NlnOperatorFas()
: hierarchy_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------*/
/* Setup of the algorithm  / operator */
void NLNSOL::NlnOperatorFas::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // create the multigrid level hierachy
  hierarchy_ = Teuchos::rcp(new NLNSOL::FAS::AMGHierarchy());
  hierarchy_->Init(Comm(), Params(), NlnProblem());
  hierarchy_->Setup();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------*/
/* Apply the preconditioner */
int NLNSOL::NlnOperatorFas::ApplyInverse(const Epetra_MultiVector& f,
    Epetra_MultiVector& x) const
{
  // Make sure that Init() and Setup() have been called
  if (not IsInit())  { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  bool converged = false;
  double fnorm2 = 1.0e+12;

  int iter = 0;
  const int itemax = Params().get<int>("max fas cycles");
  while (not converged and iter < itemax)
  {
    ++iter;

    if (Comm().MyPID() == 0)
      IO::cout << "Start V-cycle for the " << iter << ". time." << IO::endl;

    // choose type of multigrid cycle //ToDo (mayr) switch between different cycles based on params_
    VCycle(f, x, 0);

    converged = NlnProblem()->ConvergenceCheck(f, fnorm2);

    if (Comm().MyPID() == 0)
      IO::cout << IO::endl << "Finished a full VCycle run. |f| = " << fnorm2 << "." << IO::endl;
  }

  if (converged)
  {
    if (Comm().MyPID() == 0)
      IO::cout << "FAS converged after " << iter << " V-cycles." << IO::endl;
  }
  else
  {
    if (IsSolver())
      dserror("FAS failed to converge in %d V-cycles!", iter);
  }

  return 0; // ToDO (mayr) provide meaningful error code
}

/*----------------------------------------------------------------------*/
/* Do a standard V-cycle (recursive definition) */
void NLNSOL::NlnOperatorFas::VCycle(const Epetra_MultiVector& f,
    Epetra_MultiVector& x,
    const int level
    ) const
{
  int err = 0;

  Hierarchy()->CheckLevelID(level);

  if (Comm().MyPID() == 0)
    IO::cout << IO::endl << IO::endl << "WELCOME to VCycle on level " << level << IO::endl;

  // we need at least zeroed vectors, especially on the fine level
  Teuchos::RCP<Epetra_MultiVector> fbar = Teuchos::rcp(new Epetra_MultiVector(f)); // restriction of fine-level residual
  Teuchos::RCP<Epetra_MultiVector> xbar = Teuchos::rcp(new Epetra_MultiVector(x)); // restriction of fine-level solution
  Teuchos::RCP<Epetra_MultiVector> fhat = Teuchos::rcp(new Epetra_MultiVector(f.Map(), 1, true)); // coarse-grid evaluation of residual

  /* Leave 'x' untouched on the level since it is needed to approximate the
   * error after the postsmoothing. Use 'xtemp' instead for all operations on
   * this level.
   */
  Teuchos::RCP<Epetra_MultiVector> xtemp = Teuchos::rcp(new Epetra_MultiVector(x));

  // compute coarse level contributions
  if (level > 0)
  {
    // restrict fine level residual to current (coarser) level
    Hierarchy()->NlnLevel(level)->RestrictToNextCoarserLevel(fbar);

    // restrict fine level solution to current (coarser) level
    Hierarchy()->NlnLevel(level)->RestrictToNextCoarserLevel(xbar);
    Hierarchy()->NlnLevel(level)->RestrictToNextCoarserLevel(fhat);

    // set zero vectors for fhat and fbar to enable evaluation of fhat. Note: 'fhat' is zero here.
    Hierarchy()->NlnLevel(level)->NlnProblem()->SetFHatFBar(fhat, fhat);

    // evaluate coarse-grid residual
    Hierarchy()->NlnLevel(level)->NlnProblem()->Evaluate(*xbar, *fhat);

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
    dserror("Level ID %d is not a valid level ID. Level ID has to be in [0,%d]", level, Hierarchy()->NumLevels()-1);

  Teuchos::RCP<Epetra_MultiVector> fhatbar = Teuchos::rcp(new Epetra_MultiVector(*fhat));
  err = fhatbar->Update(1.0, *fbar, -1.0);
  if (err != 0) { dserror("Failed!"); }

#ifdef DEBUG
  // additional safety checks w.r.t. maps
  if (not Hierarchy()->NlnLevel(level)->DofRowMap().PointSameAs(fbar->Map())) { dserror("Maps do not match."); }
  if (not Hierarchy()->NlnLevel(level)->DofRowMap().PointSameAs(xbar->Map())) { dserror("Maps do not match."); }
  if (not Hierarchy()->NlnLevel(level)->DofRowMap().PointSameAs(fhat->Map())) { dserror("Maps do not match."); }
  if (not Hierarchy()->NlnLevel(level)->DofRowMap().PointSameAs(xtemp->Map())) { dserror("Maps do not match."); }
#endif

  // do further coarsening only in case that we are not on the coarsest level, yet.
  if (level + 1 < Hierarchy()->NumLevels())
  {
    // presmoothing
    Hierarchy()->NlnLevel(level)->DoPreSmoothing(*fhatbar, *xtemp);

    // evaluate current residual
    Teuchos::RCP<Epetra_MultiVector> fsmoothed = Teuchos::rcp(new Epetra_MultiVector(xtemp->Map(), true));
    Hierarchy()->NlnLevel(level)->NlnProblem()->Evaluate(*xtemp, *fsmoothed);

    // call VCycle on next coarser level recursively
    VCycle(*fsmoothed, *xtemp, level + 1);

    // evaluate current residual
    Hierarchy()->NlnLevel(level)->NlnProblem()->Evaluate(*xtemp, *fsmoothed);

    // postsmoothing
    Hierarchy()->NlnLevel(level)->DoPostSmoothing(*fsmoothed, *xtemp);
  }
  else // coarse level solve
  {
    if (Comm().MyPID() == 0)
    {
      IO::cout << IO::endl
               << "**************************" << IO::endl
               << "*** COARSE LEVEL SOLVE ***" << IO::endl
               << "**************************" << IO::endl
               << IO::endl;
    }

    // evaluate current residual // ToDo Do we really need to Evaluate() here?
    Teuchos::RCP<Epetra_MultiVector> fsmoothed = Teuchos::rcp(new Epetra_MultiVector(xtemp->Map(), true));
    Hierarchy()->NlnLevel(Hierarchy()->NumLevels()-1)->NlnProblem()->Evaluate(*xtemp, *fsmoothed);
    Hierarchy()->NlnLevel(Hierarchy()->NumLevels()-1)->DoCoarseLevelSolve(*fsmoothed, *xtemp);
  }

  if (level > 0)
  {
    // coarse-grid correction (=approximation of error)
    Teuchos::RCP<Epetra_MultiVector> correction = Teuchos::rcp(new Epetra_MultiVector(xtemp->Map(), 1, true));
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
  if (not x.Map().PointSameAs(f.Map())) { dserror("Map failure during recursive calls of V-cycle!"); }
#endif

  if (Comm().MyPID() == 0)
    IO::cout << "GOOD BYE from VCycle on level " << level << IO::endl;

  return;
}

/*----------------------------------------------------------------------*/
/* access the multigrid level hierarchy */
const Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy>
NLNSOL::NlnOperatorFas::Hierarchy() const
{
  if (hierarchy_.is_null())
    dserror("Hierarchy of multigrid levels 'hierarchy_' is not set, yet.");

  return hierarchy_;
}
