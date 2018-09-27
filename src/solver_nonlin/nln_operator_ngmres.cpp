/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_ngmres.cpp

\brief Nonlinear GMRES

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_Vector.h>

// standard
#include <iostream>
#include <vector>

// Teuchos
#include <Teuchos_Array.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "nln_operator_base.H"
#include "nln_operator_factory.H"
#include "nln_operator_ngmres.H"
#include "nln_problem.H"
#include "nln_utils.H"

#include "linesearch_base.H"
#include "linesearch_factory.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_sparsematrix.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorNGmres::NlnOperatorNGmres() : linesearch_(Teuchos::null), nlnprec_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorNGmres::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  // setup necessary sub-algorithms
  SetupLineSearch();
  SetupPreconditioner();

  if (getVerbLevel() > Teuchos::VERB_LOW)
  {
    *getOStream() << LabelShort() << ": Max window size has been set to " << GetMaxWindowSize()
                  << "." << std::endl;
  }

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ =
      linesearchfactory.Create(Configuration(), MyGetParameter<std::string>("line search"));

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::SetupPreconditioner()
{
  const std::string opname = MyGetParameter<std::string>("NGMRES: Nonlinear Preconditioner");

  NlnOperatorFactory nlnopfactory;
  nlnprec_ = nlnopfactory.Create(Configuration(), opname);
  nlnprec_->Init(Comm(), Configuration(), opname, NlnProblem(), BaciLinearSolver(), Nested() + 1);
  nlnprec_->Setup();

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNGmres::ApplyInverse(
    const Epetra_MultiVector& f, Epetra_MultiVector& x) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorNGmres::ApplyInverse");
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

  //----------------------------------------------------------------------------
  // prepare storage for solution iterates, corresponding residuals and norms
  //----------------------------------------------------------------------------
  std::vector<Teuchos::RCP<const Epetra_MultiVector>> sol;  // previous iterates
  std::vector<Teuchos::RCP<const Epetra_MultiVector>> res;  // previous residuals

  //----------------------------------------------------------------------------
  // some variables needed during the NGmres iteration loop
  //----------------------------------------------------------------------------
  unsigned int iter = 0;    // iteration index
  bool restart = false;     // Perform a restart?
  double fnorm2 = 1.0e+12;  // L2-norm of residual
  double fbarnorm2 = 0.0;   // L2-norm of residual after tentative step
  bool success = false;     // flag for success of acceleration step

  //----------------------------------------------------------------------------
  // criteria for acceptance of iterates and for restart decision
  // (cf. [Washio (1997)])
  //----------------------------------------------------------------------------
  // criteria
  bool criterionA = false;
  bool criterionB = false;
  bool criterionC = false;
  bool criterionD = false;
  bool criterionCprev = false;
  bool criterionDprev = false;

  // constants
  const double gammaA = 2.0;
  const double epsB = 0.1;
  const double deltaB = 0.9;
  const double gammaC = std::max(2.0, gammaA);
  //----------------------------------------------------------------------------

  if (getVerbLevel() > Teuchos::VERB_LOW)
  {
    *getOStream() << "Begin with Krylov acceleration, now. "
                  << "Starting iteration count is " << iter << "." << std::endl;
  }

  //----------------------------------------------------------------------------
  // global solution and residual vectors
  //----------------------------------------------------------------------------
  // current iterate
  Teuchos::RCP<Epetra_MultiVector> xbar = Teuchos::rcp(new Epetra_MultiVector(x));

  // residual at current iterate
  Teuchos::RCP<Epetra_MultiVector> fbar = Teuchos::rcp(new Epetra_MultiVector(f.Map(), true));
  NlnProblem()->ComputeF(*xbar, *fbar);

  // accelerated iterate
  Teuchos::RCP<Epetra_MultiVector> xhat = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));

  // residual at accelerated iterate
  Teuchos::RCP<Epetra_MultiVector> fhat = Teuchos::rcp(new Epetra_MultiVector(f.Map(), true));
  //----------------------------------------------------------------------------

  bool converged = NlnProblem()->ConvergenceCheck(*fbar, fnorm2);

  //----------------------------------------------------------------------------
  // setup stagnation detection mechanism
  //----------------------------------------------------------------------------
  Teuchos::RCP<NLNSOL::UTILS::StagnationDetection> stagdetect =
      Teuchos::rcp(new NLNSOL::UTILS::StagnationDetection());
  stagdetect->Init(Configuration(),
      MyGetParameter<std::string>("nonlinear operator: stagnation detection"), fnorm2);
  //----------------------------------------------------------------------------

  // print initial state of convergence
  PrintIterSummary(iter, fnorm2);

  //----------------------------------------------------------------------------
  // outer iteration loop
  //----------------------------------------------------------------------------
  while (ContinueIterations(iter, converged))
  {
    // reset the restart flag
    restart = false;

    //--------------------------------------------------------------------------
    // window iteration loop
    //--------------------------------------------------------------------------
    while (not restart)
    {
      ++iter;

      // add most recent iterate and residual to the history
      AddToWindow(xbar, sol);
      AddToWindow(fbar, res);

      //------------------------------------------------------------------------
      // Step 1: Generate a new tentative iterate
      //------------------------------------------------------------------------
      {
        // compute a tentative iterate by applying the preconditioner once
        ComputeTentativeIterate(*fbar, *xbar);

        NlnProblem()->ComputeF(*xbar, *fbar);
        converged = NlnProblem()->ConvergenceCheck(*fbar, fbarnorm2);
        if (converged) break;
      }
      //------------------------------------------------------------------------

      //------------------------------------------------------------------------
      // Step 2: Generate accelerated iterate
      //------------------------------------------------------------------------
      {
        success = ComputeAcceleratedIterate(xbar, fbar, sol, res, xhat);
        if (not success)
        {
          restart = true;

          if (getVerbLevel() > Teuchos::VERB_LOW)
            *getOStream() << "*** Acceleration failed in iteration    " << iter << std::endl;
        }

        // evaluate the residual based on the accelerated solution
        NlnProblem()->ComputeF(*xhat, *fhat);
      }
      //------------------------------------------------------------------------

      //------------------------------------------------------------------------
      // Evaluate criteria A, B, C, and D as in [Washio (1997)]
      //------------------------------------------------------------------------
      {
        // update
        criterionCprev = criterionC;
        criterionDprev = criterionD;
        criterionC = false;
        criterionD = false;

        //----------------------------------------------------------------------
        // Compute norms
        //----------------------------------------------------------------------
        double fnormhat = 0.0;
        NlnProblem()->ConvergenceCheck(*fhat, fnormhat);

        double fnormmin = fbarnorm2;
        for (unsigned int i = 0; i < res.size(); ++i)
        {
          double tmpnorm = 1.0e+12;
          NlnProblem()->ConvergenceCheck(*(res[i]), tmpnorm);
          fnormmin = std::min(fnormmin, tmpnorm);
        }

        Teuchos::RCP<Epetra_MultiVector> vec =
            Teuchos::rcp(new Epetra_MultiVector(xhat->Map(), true));
        double xnormmin = 1.0e+12;
        for (unsigned int i = 0; i < sol.size(); ++i)
        {
          double tmpnorm = 1.0e+12;
          vec->Update(1.0, *xhat, -1.0, *(sol[i]), 0.0);
          NlnProblem()->ConvergenceCheck(*vec, tmpnorm);
          xnormmin = std::min(xnormmin, tmpnorm);
        }

        double xnormleft = 0.0;
        vec->Update(1.0, *xhat, -1.0, *xbar, 0.0);
        NlnProblem()->ConvergenceCheck(*vec, xnormleft);

        //----------------------------------------------------------------------
        // Check all criteria
        //----------------------------------------------------------------------
        // Criterion A
        if (fnormhat < gammaA * fnormmin) criterionA = true;

        // criterion B
        if (epsB * xnormleft < xnormmin or fnormhat < deltaB * fnormmin) criterionB = true;

        // Criterion C
        if (fnormhat > gammaC * fnormmin) criterionC = true;

        // Criterion D
        if (epsB * xnormleft >= xnormmin and fnormhat >= deltaB * fnormmin) criterionD = true;

        //----------------------------------------------------------------------
        // Restart necessary?
        //----------------------------------------------------------------------
        if (not success or (criterionC and criterionCprev) or (criterionD and criterionDprev))
        {
          restart = true;
        }

        //----------------------------------------------------------------------
        // Accept current iterate?
        //----------------------------------------------------------------------
        if (criterionA and criterionB and not restart)
        {
          xbar->Update(1.0, *xhat, 0.0);
          NlnProblem()->ComputeF(*xbar, *fbar);
        }
      }
      //------------------------------------------------------------------------

      //------------------------------------------------------------------------
      // Step 3: Perform line search as in [Sterck (2012)]
      //------------------------------------------------------------------------
      //      {
      //        if (not success
      //            or (criterionC and criterionD and criterionCprev and criterionDprev))
      //        {
      //          restart = true;
      //        }
      //        else // accepted iterate and, thus, perform a line search step
      //        {
      //          PerformLineSearchStep(xbar, xhat, fbar);
      //          NlnProblem()->ComputeF(*xbar, *fbar);
      //        }
      //      }
      //------------------------------------------------------------------------

      if (not(sol.size() <=
              GetMaxWindowSize()))  // window size reached maximum ToDo (mayr) needs to be <=
      {
        // erase the oldest solution and residual vectors
        sol.erase(sol.begin());
        res.erase(res.begin());
      }

      // perform the restart
      if (restart)
      {
        if (getVerbLevel() > Teuchos::VERB_LOW)
        {
          *getOStream() << LabelShort() << ": Perform restart in iteration    " << iter
                        << std::endl;
        }

        // clear history
        sol.clear();
        res.clear();

        // reset restart criteria
        criterionC = false;
        criterionD = false;
        criterionCprev = false;
        criterionDprev = false;
      }

      if (getVerbLevel() > Teuchos::VERB_LOW)
        *getOStream() << LabelShort() << ": Current window size:   " << sol.size() << std::endl;

      converged = NlnProblem()->ConvergenceCheck(*fbar, fnorm2);

      // check for stagnation
      if (stagdetect->Check(fnorm2))
      {
        *getOStream() << LabelShort() << " detected stagnation. Rebuild preconditioner."
                      << std::endl;
        nlnprec_->RebuildPrec();
      }

      // print stuff
      PrintIterSummary(iter, fnorm2);

      if (iter > (unsigned int)GetMaxIter() or converged) break;
    }  // end of while for window
  }    // end of while for outer iteration loop

  // update solution that is passed back to the calling algorithm
  err = x.Update(1.0, *xbar, 0.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  // final convergence check and output
  NlnProblem()->ComputeF(x, *fbar);
  NlnProblem()->ConvergenceCheck(*fbar, fnorm2);
  PrintIterSummary(iter, fnorm2);

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
  SetOutParameterStagnation(stagdetect->StatusParams());

  // return error code
  return errorcode;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNGmres::ComputeTentativeIterate(
    const Epetra_MultiVector& fbar, Epetra_MultiVector& xbar) const
{
  dsassert(xbar.Map().PointSameAs(fbar.Map()), "Maps do not match.");

  int errorcode = nlnprec_->ApplyInverse(fbar, xbar);

  // Do we need to rebuild the preconditioner?
  if (nlnprec_->GetOutParams()->get<NLNSOL::UTILS::OperatorStatus>("Error Code") ==
      NLNSOL::UTILS::opstatus_stagnation)
  {
    // create formatted output stream
    Teuchos::RCP<Teuchos::FancyOStream> out =
        Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    Teuchos::OSTab tab(out, Indentation());

    *getOStream() << "Preconditioner of " << LabelShort() << " needs to be rebuild." << std::endl;

    nlnprec_->RebuildPrec();
  }

  return errorcode;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnOperatorNGmres::ComputeAcceleratedIterate(
    const Teuchos::RCP<const Epetra_MultiVector> xbar,
    const Teuchos::RCP<const Epetra_MultiVector> fbar,
    const std::vector<Teuchos::RCP<const Epetra_MultiVector>>& sol,
    const std::vector<Teuchos::RCP<const Epetra_MultiVector>>& res,
    Teuchos::RCP<Epetra_MultiVector> xhat) const
{
  int err = 0;

  bool success = true;  // suppose success of acceleration step

  if (sol.size() != res.size() or sol.size() < 1 or res.size() < 1)
    dserror("Size of history vectors 'sol' or 'res' is not admissible.");

  // create quantities needed for normal equations of least squares problem
  LINALG::SerialDenseVector alpha;  // weights for linear combination of Krylov vectors
  LINALG::SerialDenseVector Ptg;    // right hand side vector
  LINALG::SerialDenseMatrix PtP;    // "system matrix"

  // set suitable sizes for least squares problem and initialize to zero
  {
    const int wsize = sol.size();  // current window size
    alpha.Size(wsize);
    Ptg.Size(wsize);
    PtP.Shape(wsize, wsize);
  }

  // ---------------------------------------------------------------------------
  // construct least squares problem (cf. [Sterck (2012), eq. (2.7)])
  // ---------------------------------------------------------------------------
  // matrix P
  std::vector<Teuchos::RCP<Epetra_MultiVector>> P;
  P.clear();
  {
    for (std::vector<Teuchos::RCP<const Epetra_MultiVector>>::const_iterator it = res.begin();
         it < res.end(); ++it)
    {
      // allocate new vector as column of matrix P
      Teuchos::RCP<Epetra_MultiVector> p = Teuchos::rcp(new Epetra_MultiVector(fbar->Map(), true));

      /* difference between most recent residual and the one from the 'it'
       * previous iteration */
      err = p->Update(1.0, *fbar, -1.0, *(*it), 0.0);
      if (err != 0)
      {
        dserror("Update failed.");
      }

      // insert this as last column into matrix 'P'
      P.push_back(p);
    }
  }

  // do matrix-matrix product P^T * P and right hand side product P^T * g
  double maxdiagentry =
      0.0;  // max value on main diagonal (needed to fix possible singularity of 'P^T * P')
  {
    double rhsentry = 0.0;     // entry to right hand side at row (i)
    double matrixentry = 0.0;  // entry to matrix at location (i,j)
    for (unsigned int i = 0; i < P.size(); ++i)
    {
      // right hand side product -P^T * g
      P[i]->Dot(*fbar, &rhsentry);
      Ptg[i] = -rhsentry;

      // matrix-matrix product P^T * P
      for (unsigned int j = 0; j <= i; ++j)
      {
        // compute single entry of target matrix PtP
        P[i]->Dot(*(P[j]), &matrixentry);

        // put matrix entry at its positions (symmetric matrix)
        if (i != j)  // off-diagonal entries (symmetry!)
        {
          PtP(i, j) = matrixentry;
          PtP(j, i) = matrixentry;
        }
        else  // main diagonal entries
        {
          PtP(i, i) = matrixentry;

          // determine max main diagonal entry
          maxdiagentry = std::max(matrixentry, maxdiagentry);
        }
      }
    }
  }

  /* Fix diagonal of 'PtP' to ensure non-singularity
   *
   * Note: There is no theoretical reason why we need this. It's more a
   * technical fix of numerical problems that might occur sometimes. However,
   * [Washio (1997), Lemma 2.2] confirms that this modification produces only
   * negligible errors in the coefficients.
   */
  const double delta = 1.0e-12 * maxdiagentry;
  for (int i = 0; i < PtP.RowDim(); ++i)
  {
    PtP(i, i) += delta;
  }

  {
    /* Solve normal equations with direct solver since it is just a very small
     * system (size = window size) */

    // create solver
    Epetra_SerialDenseSolver solver;

    // set the linear system
    solver.SetMatrix(PtP);
    solver.SetVectors(alpha, Ptg);

    // configure solver strategy
    solver.FactorWithEquilibration(true);
    solver.SolveToRefinedSolution(true);

    if (solver.Factor() or solver.Solve()) success = false;
  }

  // ---------------------------------------------------------------------------

  // sum up last iterates to new solution, i.e. do the linear combination
  err = xhat->Update(1.0, *xbar, 0.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }
  if (success)
  {
    for (int i = 0; i < alpha.Length(); ++i)
    {
      err = xhat->Update(alpha[i], *xbar, -alpha[i], *(sol[i]), 1.0);
      if (err != 0)
      {
        dserror("Update failed.");
      }
    }
  }

  return success;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::PerformLineSearchStep(Teuchos::RCP<Epetra_MultiVector> xbar,
    Teuchos::RCP<const Epetra_MultiVector> xhat, Teuchos::RCP<Epetra_MultiVector> fbar) const
{
  // some initializations
  double resnormold = 0.0;  // residual norm at current iterate
  double steplength = 1.0;  // line search parameter
  bool suffdecr = false;    // flag for sufficient decrease of line search

  // build vector of search direction
  Teuchos::RCP<Epetra_MultiVector> inc = Teuchos::rcp(new Epetra_MultiVector(*xhat));
  inc->Update(-1.0, *xbar, 1.0);

  // line search
  NlnProblem()->ConvergenceCheck(*fbar, resnormold);
  ComputeStepLength(*xbar, *fbar, *inc, resnormold, steplength, suffdecr);

  // update solution using the line search parameter
  xbar->Update(steplength, *xhat, 1.0 - steplength);

  return;
}

/*----------------------------------------------------------------------------*/
unsigned int NLNSOL::NlnOperatorNGmres::GetMaxWindowSize() const
{
  return MyGetParameter<int>("NGMRES: max window size");
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::ComputeStepLength(const Epetra_MultiVector& x,
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
void NLNSOL::NlnOperatorNGmres::AddToWindow(Teuchos::RCP<const Epetra_MultiVector> vec,
    std::vector<Teuchos::RCP<const Epetra_MultiVector>>& history) const
{
  Teuchos::RCP<const Epetra_MultiVector> veccopy = Teuchos::rcp(new Epetra_MultiVector(*vec));
  history.push_back(veccopy);

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::RebuildPrec()
{
  nlnprec_->RebuildPrec();

  return;
}
