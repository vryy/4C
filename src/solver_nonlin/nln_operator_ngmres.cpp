/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_ngmres.cpp

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
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_Vector.h>

// standard
#include <iostream>
#include <vector>

// Teuchos
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_base.H"
#include "nln_operator_factory.H"
#include "nln_operator_ngmres.H"
#include "nln_problem.H"

#include "linesearch_base.H"
#include "linesearch_factory.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorNGmres::NlnOperatorNGmres()
: linsolver_(Teuchos::null),
  linesearch_(Teuchos::null),
  nlnprec_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // setup necessary sub-algorithms
  SetupLinearSolver();
  SetupLineSearch();
  SetupPreconditioner();

  if (Comm().MyPID() == 0)
    IO::cout << "Max window size has been set to " << GetMaxWindowSize()
             << "." << IO::endl;

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::SetupLinearSolver()
{
  // get the solver number used for structural problems
  const int linsolvernumber = Params().get<int>("NGMRES: Linear Solver");

  // check if the solver ID is valid
  if (linsolvernumber == (-1))
    dserror("No valid linear solver defined!");

  linsolver_ = Teuchos::rcp(
      new LINALG::Solver(
          DRT::Problem::Instance()->SolverParams(linsolvernumber), Comm(),
          DRT::Problem::Instance()->ErrorFile()->Handle()));

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ = linesearchfactory.Create(
      Params().sublist("NGMRES: Line Search"));

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNGmres::SetupPreconditioner()
{
  const Teuchos::ParameterList& precparams =
      Params().sublist("NGMRES: Nonlinear Preconditioner");

  NlnOperatorFactory nlnopfactory;
  nlnprec_ = nlnopfactory.Create(precparams);
  nlnprec_->Init(Comm(), precparams, NlnProblem());
  nlnprec_->Setup();

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNGmres::ApplyInverse(const Epetra_MultiVector& f,
    Epetra_MultiVector& x
    ) const
{
  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // prepare storage for solution iterates and corresponding residuals
  std::vector<Teuchos::RCP<Epetra_MultiVector> > sol;
  std::vector<Teuchos::RCP<Epetra_MultiVector> > res;

  // some variables needed during the NGmres iteration loop
  unsigned int iter = 0; // iteration index
  double steplength = 1.0; // line search parameter
  double resnormold = 1.0e+12; // residual L2 norm
  bool restart = false; // Perform a restart?
  double dotproduct = 0.0; // dot product of residual with search direction
  double fnorm2 = 1.0e+12; // L2-norm of residual

  if (Comm().MyPID() == 0)
  {
    IO::cout << "Begin with Krylov acceleration, now. "
             << "Starting iteration count is " << iter << "." << IO::endl;
  }

  // get local copies of solution and residual
  Teuchos::RCP<Epetra_MultiVector> xbar =
      Teuchos::rcp(new Epetra_MultiVector(x));
  Teuchos::RCP<Epetra_MultiVector> fbar =
      Teuchos::rcp(new Epetra_MultiVector(f.Map(), true));
  NlnProblem()->ComputeF(*xbar, *fbar);

  bool converged = NlnProblem()->ConvergenceCheck(*fbar, fnorm2);

  // print stuff
  PrintIterSummary(iter, fnorm2);

  // outer iteration loop
  while (ContinueIterations(iter, converged))
  {
    // reset the restart flag
    restart = false;

    // window iteration loop
    while (not restart)
    {
      ++iter;

      // add most recent iterate and residual to the history
      {
        Teuchos::RCP<Epetra_MultiVector> xcopy =
            Teuchos::rcp(new Epetra_MultiVector(*xbar));
        Teuchos::RCP<Epetra_MultiVector> fcopy =
            Teuchos::rcp(new Epetra_MultiVector(*fbar));
        sol.push_back(xcopy);
        res.push_back(fcopy);
      }

      /* -------------------------------------------------------------------- */
      // Step 1: Generate a new tentative iterate
      /* -------------------------------------------------------------------- */
      err = xbar->Update(1.0, *sol.back(), 0.0);
      if (err != 0) { dserror("Update failed."); }
      err = fbar->Update(1.0, *res.back(), 0.0);
      if (err != 0) { dserror("Update failed."); }

      // compute a tentative iterate by applying the preconditioner once
      ComputeTentativeIterate(*fbar, *xbar);

      // evaluate the residual based on the new tentative solution
      NlnProblem()->ComputeF(*xbar, *fbar);
      /* -------------------------------------------------------------------- */

      /* -------------------------------------------------------------------- */
      // Step 2: Generate accelerated iterate
      /* -------------------------------------------------------------------- */
      Teuchos::RCP<Epetra_MultiVector> xhat =
          ComputeAcceleratedIterate(xbar, fbar, sol, res);

      /* -------------------------------------------------------------------- */

      /* -------------------------------------------------------------------- */
      // Step 3: Perform line search
      /* -------------------------------------------------------------------- */
      // the full step increment
      Teuchos::RCP<Epetra_MultiVector> inc =
          Teuchos::rcp(new Epetra_MultiVector(*xhat));
      err = inc->Update(-1.0, *xbar, 1.0);
      if (err != 0) { dserror("Update failed."); }

      // Is 'inc' a descent direction? Otherwise restart the algorithm.
      inc->Dot(*fbar, &dotproduct);
      if (dotproduct >= 0.0) // ascent direction --> restart required
      {
        restart = true;
        if (Comm().MyPID() == 0)
          IO::cout << "Perform restart in iteration " << iter << "."
                   << IO::endl;
      }
      else // descent direction
      {
        // line search
        NlnProblem()->ConvergenceCheck(*fbar, resnormold);
        steplength = ComputeStepLength(*xbar, *inc, resnormold);

        /* update solution using the line search parameter
         * (called 'xstar' in [Sterck2012a]) */
        err = xbar->Update(steplength, *xhat, 1.0-steplength);
        if (err != 0) { dserror("Update failed."); }

        // evaluate the residual based on the new solution
        NlnProblem()->ComputeF(*xbar, *fbar);
      }
      /* -------------------------------------------------------------------- */

      if (not (sol.size() < GetMaxWindowSize())) // window size reached maximum
      {
        // erase the oldest solution and residual vectors
        sol.erase(sol.begin());
        res.erase(res.begin());
      }

      if (restart)
      {
        sol.clear();
        res.clear();
      }

      if (Comm().MyPID() == 0)
        IO::cout << "Current window size:   " << sol.size() << IO::endl;

      converged = NlnProblem()->ConvergenceCheck(*fbar, fnorm2);

      // print stuff
      PrintIterSummary(iter, fnorm2);

      if (iter > (unsigned int)GetMaxIter() or converged)
        break;
    } // end of while for window
  } // end of while for outer iteration loop

  // update solution that is passed back to the calling algorithm
  err = x.Update(1.0, *xbar, 0.0);
  if (err != 0) { dserror("Update failed."); }

  // return error code
  return (not CheckSuccessfulConvergence(iter, converged));
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNGmres::ComputeTentativeIterate(
    const Epetra_MultiVector& fbar, Epetra_MultiVector& xbar) const
{
  if (not xbar.Map().PointSameAs(fbar.Map()))
    dserror("Maps do not match.");

  return nlnprec_->ApplyInverse(fbar, xbar);
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector>
NLNSOL::NlnOperatorNGmres::ComputeAcceleratedIterate(
    const Teuchos::RCP<const Epetra_MultiVector>& xbar,
    const Teuchos::RCP<const Epetra_MultiVector>& fbar,
    const std::vector<Teuchos::RCP<Epetra_MultiVector> >& sol,
    const std::vector<Teuchos::RCP<Epetra_MultiVector> >& res
    ) const
{
  int err = 0;

  // create quantities needed for least squares problem
  LINALG::SerialDenseVector alpha; // weights for linear combination
  LINALG::SerialDenseVector Ptg; // right hand side vector
  LINALG::SerialDenseMatrix PtP; // "system matrix"

  // set suitable sizes for least squares problem and initialize to zero
  const int wsize = sol.size(); // current window size
  alpha.Size(wsize);
  Ptg.Size(wsize);
  PtP.Shape(wsize, wsize);

  // construct least squares problem
  Teuchos::Array<Teuchos::RCP<Epetra_MultiVector> > P;
  P.clear();
  {
    Teuchos::RCP<Epetra_MultiVector> p =
        Teuchos::rcp(new Epetra_MultiVector(fbar->Map(), true));
    for (std::vector<Teuchos::RCP<Epetra_MultiVector> >::const_iterator it = res.begin(); it < res.end(); ++it)
    {
      /* difference between most recent residual and the one from the 'it'
       * previous iteration */
      err = p->Update(1.0, *fbar, -1.0, *(*it), 0.0);
      if (err != 0) { dserror("Update failed."); }

      // insert this as last column into matrix 'P'
      P.push_back(p);
    }
  }

  // do matrix-matrix product P^T * P and right hand side product P^T * g
  double maxdiagentry = 0.0; // needed to fix possible singularity of 'P^T * P'
  {
    double rhsentry = 0.0;
    double matrixentry = 0.0;
    for (unsigned int i = 0; i < P.size(); ++i)
    {
      // right hand side product P^T * g
      err = P[i]->Dot(*fbar, &rhsentry);
      if (err != 0) { dserror("Failed."); }
      Ptg[i] = -1.0 * rhsentry;

      // matrix-matrix product P^T * P
      for (unsigned int j = 0; j <= i; ++j)
      {
        // compute single entry of target matrix PtP
        err = P[i]->Dot(*P[j], &matrixentry);
        if (err != 0) { dserror("Failed."); }

        // put matrix entry at its positions (symmetric matrix)
        if (i != j) // off-diagonal entries (symmetry!)
        {
          PtP[i][j] = matrixentry;
          PtP[j][i] = matrixentry;
        }
        else // main diagonal entries
        {
          PtP[i][i] = matrixentry;

          // determine max main diagonal entry
          maxdiagentry = std::max(matrixentry, maxdiagentry);
        }
      }
    }
  }

  // fix diagonal of 'PtP' to ensure non-singularity
  /* Literature proposed different options how to ensure the non-singularity
   * of the matrix PtP. Currently, we have too less experience to prefer one
   * or the other.
   *
   * ToDo (mayr) Drop one of them or introduce control via input flag
   */
  for (int i = 0; i < PtP.RowDim(); ++i)
  {
    PtP[i][i] += 1.0e-8 * maxdiagentry; // [Sterck2012a]
//    PtP[i][i] += 1.0e-6; // [Washio1997a]
  }

  /* Solve normal equations with direct solver since it is just a very small
   * system (size = window size) */
  {
    // create solver
    Epetra_SerialDenseSolver solver;

    // set the linear system
    solver.SetMatrix(PtP);
    solver.SetVectors(alpha, Ptg);

    // configure solver strategy
    solver.FactorWithEquilibration(true);
    solver.SolveToRefinedSolution(true);

    // solve
    int err1 = solver.Factor();
    int err2 = solver.Solve();

    // check for errors
    if ( err1 != 0 or err2 != 0 )
      dserror("Something went wrong! Factor() and Solve() gave error codes %d "
          "and %d", err1, err2);
  }

  // sum up last iterates to new solution, i.e. do the linear combination
  Teuchos::RCP<Epetra_MultiVector> xhat =
      Teuchos::rcp(new Epetra_MultiVector(*xbar));
  if (xhat.is_null()) { dserror("Allocation failed."); }
  for (int i = 0; i < alpha.Length(); ++i)
  {
    err = xhat->Update(alpha[i], *xbar, -alpha[i], *(sol[i]), 1.0);
    if (err != 0) { dserror("Update failed."); }
  }

  return xhat;
}

/*----------------------------------------------------------------------------*/
const unsigned int NLNSOL::NlnOperatorNGmres::GetMaxWindowSize() const
{
  return Params().get<int>("NGMRES: Max Window Size");
}

/*----------------------------------------------------------------------------*/
const double NLNSOL::NlnOperatorNGmres::ComputeStepLength(
    const Epetra_MultiVector& x, const Epetra_MultiVector& inc,
    double fnorm2) const
{
  linesearch_->Init(NlnProblem(), Params().sublist("NGMRES: Line Search"), x,
      inc, fnorm2);
  linesearch_->Setup();
  return linesearch_->ComputeLSParam();
}
