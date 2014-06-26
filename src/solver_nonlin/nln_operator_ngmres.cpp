/*----------------------------------------------------------------------*/
/*!
\file nln_operator_ngmres.cpp

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
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_Vector.h>

// standard
#include <iostream>
#include <vector>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_factory.H"
#include "nln_operator_ngmres.H"
#include "nln_operator.H"
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

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperatorNGmres::NlnOperatorNGmres()
: linsolver_(Teuchos::null),
  linesearch_(Teuchos::null),
  nlnprec_(Teuchos::null),
  itermax_(0),
  windowsize_(0)
{
  return;
}

/*----------------------------------------------------------------------*/
/* Setup of the algorithm */
void NLNSOL::NlnOperatorNGmres::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // setup necessary sub-algorithms
  SetupLinearSolver();
  SetupLineSearch();
  SetupPreconditioner();

  itermax_ = Params().get<int>("NGMRES: Max Iter");
  windowsize_ = Params().get<int>("NGMRES: Max Window Size");

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------*/
/* Setup of the linear solver */
void NLNSOL::NlnOperatorNGmres::SetupLinearSolver()
{
  // get the solver number used for structural problems
  const int linsolvernumber = Params().get<int>("NGMRES: Linear Solver");

  // check if the solver ID is valid
  if (linsolvernumber == (-1))
    dserror("No valid linear solver defined!");

  linsolver_ = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                               Comm(),
                                               DRT::Problem::Instance()->ErrorFile()->Handle()));

  return;
}

/*----------------------------------------------------------------------*/
/* Setup of the line search */
void NLNSOL::NlnOperatorNGmres::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ = linesearchfactory.Create(Params().sublist("NGMRES: Line Search"));

  return;
}

/*----------------------------------------------------------------------*/
/* Setup of the preconditioner */
void NLNSOL::NlnOperatorNGmres::SetupPreconditioner()
{
  const Teuchos::ParameterList& precparams = Params().sublist("NGMRES: Nonlinear Preconditioner");

  NlnOperatorFactory nlnopfactory;
  nlnprec_ = nlnopfactory.Create(precparams);
  nlnprec_->Init(Comm(), precparams, NlnProblem());
  nlnprec_->Setup();

  return;
}

/*----------------------------------------------------------------------*/
/* ApplyInverse of nonlinear Krylov solver */
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

  // save initial solution and initial residual
  {
    Teuchos::RCP<Epetra_MultiVector> xcopy = Teuchos::rcp(new Epetra_MultiVector(x));
    Teuchos::RCP<Epetra_MultiVector> fcopy = Teuchos::rcp(new Epetra_MultiVector(f.Map(), true));
    NlnProblem()->Evaluate(*xcopy, *fcopy);
    sol.push_back(xcopy);
    res.push_back(fcopy);
  }

  // Krylov loop iteration index and its upper and lower bounds
  int iter = 0;

  double steplength = 1.0; // line search parameter
  double resnormold = 1.0e+12; // residual L2 norm

  if (Comm().MyPID() == 0)
  {
    IO::cout << "Begin with Krylov acceleration, now. "
             << "Starting iteration count is " << iter << "." << IO::endl;
  }

  // declare some vectors
  Teuchos::RCP<Epetra_MultiVector> xbar = Teuchos::rcp(new Epetra_MultiVector(x));
  Teuchos::RCP<Epetra_MultiVector> fbar = Teuchos::rcp(new Epetra_MultiVector(f.Map(), true));
  NlnProblem()->Evaluate(*xbar, *fbar);

  bool converged =  NlnProblem()->ConvergenceCheck(*fbar);

  if (not xbar->Map().PointSameAs(fbar->Map()))
    dserror("Maps do not match.");

  while ((not converged) and (iter <= itermax_))
  {
    // start a new window
    {
      // store initial solution and residual
      Teuchos::RCP<Epetra_MultiVector> xcopy = Teuchos::rcp(new Epetra_MultiVector(*xbar));
      Teuchos::RCP<Epetra_MultiVector> fcopy = Teuchos::rcp(new Epetra_MultiVector(*fbar));
      sol.push_back(xcopy);
      res.push_back(fcopy);
    }

    bool restart = false; // Perform a restart?
    int wsize = 1; // current window size

    while (not restart)
    {
      ++iter;

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
      NlnProblem()->Evaluate(*xbar, *fbar);
      /* -------------------------------------------------------------------- */

      /* -------------------------------------------------------------------- */
      // Step 2: Solve least squares problem and do the linear combination
      /* -------------------------------------------------------------------- */
      Teuchos::RCP<Epetra_MultiVector> xhat = ComputeAcceleratedIterate(xbar, fbar, sol, res, wsize);

      /* -------------------------------------------------------------------- */

      /* -------------------------------------------------------------------- */
      // Step 3: Perform line search
      /* -------------------------------------------------------------------- */
      // the full step increment
      Teuchos::RCP<Epetra_MultiVector> inc = Teuchos::rcp(new Epetra_MultiVector(*xhat));
      err = inc->Update(-1.0, *xbar, 1.0);
      if (err != 0) { dserror("Update failed."); }

      // Is 'inc' a descent direction? Otherwise restart the algorithm.
      double dotproduct = 0.0;
      inc->Dot(*fbar, &dotproduct);
      if (dotproduct >= 0.0)
      {
        restart = true;
        if (Comm().MyPID() == 0)
          IO::cout << "Perform restart in iteration " << iter << "." << IO::endl;
      }
      else
      {
        // line search
        NlnProblem()->ConvergenceCheck(*fbar, resnormold);

        linesearch_->Init(NlnProblem(), Params().sublist("NGMRES: Line Search"), *xbar, *inc, resnormold);
        linesearch_->Setup();
        steplength = linesearch_->ComputeLSParam();

        // update solution using the line search parameter
        Teuchos::RCP<Epetra_MultiVector> xstar = Teuchos::rcp(new Epetra_MultiVector(*xhat));
        err = xstar->Update(1.0-steplength, *xbar, steplength);
        if (err != 0) { dserror("Update failed."); }

        // evaluate the residual based on the new tentative solution
        NlnProblem()->Evaluate(*xstar, *fbar);

        {
          Teuchos::RCP<Epetra_MultiVector> xcopy = Teuchos::rcp(new Epetra_MultiVector(*xstar));
          Teuchos::RCP<Epetra_MultiVector> fcopy = Teuchos::rcp(new Epetra_MultiVector(*fbar));
          err = fcopy->Norm2(&resnormold);
          if (err != 0) { dserror("Failed!"); }
          resnormold /= sqrt(fcopy->GlobalLength());

          sol.push_back(xcopy);
          res.push_back(fcopy);
        }

        x.Update(1.0, *xstar, 0.0);
      }
      /* -------------------------------------------------------------------- */

      if (wsize <= windowsize_) // window is still growing
      {
        ++wsize;
      }
      else // window size is at maximum, so we delete the oldest values
      {
        sol.erase(sol.begin());
        res.erase(res.begin());
      }

      if (restart)
      {
        sol.clear();
        res.clear();
        wsize = 1;
      }

      if (Comm().MyPID() == 0)
        IO::cout << "Current window size:   " << wsize << IO::endl;

      double fnorm2 = 1.0e+12;
      converged = NlnProblem()->ConvergenceCheck(*fbar, fnorm2);

      // print stuff
      PrintIterSummary(iter, fnorm2);

      if ( iter > itermax_ or converged)
        break;
    } // end of while for window
  } // end of while for outer iteration loop

  if (IsSolver() and iter > itermax_ and not converged)
  {
    dserror("Nonlinear Krylov acceleration did not converge in %i iterations", itermax_);
  }

  // return error code
  return 0; // ToDo (mayr) return meaningful error code
}

/*----------------------------------------------------------------------*/
/* Compute a new tentative iterate */
int NLNSOL::NlnOperatorNGmres::ComputeTentativeIterate(const Epetra_MultiVector& fbar,
    Epetra_MultiVector& xbar
    ) const
{
  return nlnprec_->ApplyInverse(fbar, xbar);
}

/*----------------------------------------------------------------------*/
/* Solve the least squares problem */
Teuchos::RCP<Epetra_MultiVector> NLNSOL::NlnOperatorNGmres::ComputeAcceleratedIterate
(
  const Teuchos::RCP<Epetra_MultiVector>& xbar,
  const Teuchos::RCP<Epetra_MultiVector>& fbar,
  const std::vector<Teuchos::RCP<Epetra_MultiVector> >& sol,
  const std::vector<Teuchos::RCP<Epetra_MultiVector> >& res,
  const int wsize
) const
{
  int err = 0;

  // create quantities needed for least squares problem
  LINALG::SerialDenseVector alpha;
  LINALG::SerialDenseVector Ptg;
  LINALG::SerialDenseMatrix PtP;

  // set suitable sizes for least squares problem and initialize to zero
  alpha.Size(wsize);
  Ptg.Size(wsize);
  PtP.Shape(wsize, wsize);

  // construct least squares problem
  std::vector<Teuchos::RCP<Epetra_MultiVector> > P;
  P.clear();
  for (int colidx = 0; colidx < wsize; ++colidx)
  {
    // difference between most recent residual and the one from iteration 'colidx'
    Teuchos::RCP<Epetra_MultiVector> p = Teuchos::rcp(new Epetra_MultiVector(*fbar));
    err = p->Update(-1.0, *(res[colidx]), 1.0);
    if (err != 0) { dserror("Update failed."); }

    // put this into matrix 'P' into column 'colidx'
    P.push_back(p);
  }

  double maxdiagentry = 0.0; // needed to fix possible singularity of 'P^y * P'

  // do matrix-matrix product P^T * P and right hand side product P^T * g
  {
    double rhsentry = 0.0;
    double matrixentry = 0.0;
    for (int i = 0; i < wsize; ++i)
    {
      // right hand side
      err = P[i]->Dot(*fbar, &rhsentry);
      if (err != 0) { dserror("Failed."); }
      Ptg[i] = -1.0 * rhsentry;

      // matrix-matrix product
      for (int j = 0; j <= i; ++j)
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
        else        // main diagonal entries
        {
          PtP[i][i] = matrixentry;
          if (matrixentry > maxdiagentry) // determine max main diagonal entry
            maxdiagentry = matrixentry;
        }
      }
    }
  }

  // fix diagonal of 'PtP' to ensure non-singularity
  for (int i = 0; i < wsize; ++i)
  {
    PtP[i][i] += 1.0e-8 * maxdiagentry; // [Sterck2012a]
//    PtP[i][i] += 1.0e-6; // [Washio1997a]
  }

  // Solve normal equations with direct solver since it is just a very small system
  {
    Epetra_SerialDenseSolver solver;
    solver.SetMatrix(PtP);
    solver.SetVectors(alpha, Ptg);
    solver.FactorWithEquilibration(true);
    int err1 = solver.Factor();
    int err2 = solver.Solve();
    if ( err1 != 0 or err2 != 0 )
      dserror("Something went wrong! Factor() and Solve() gave error codes %d and %d", err1, err2);
  }

  // sum up last iterates to new solution, i.e. do the linear combination
  Teuchos::RCP<Epetra_MultiVector> xhat = Teuchos::rcp(new Epetra_MultiVector(*xbar));
  if (xhat.is_null()) { dserror("Allocation failed."); }
  for (int i = 0; i < alpha.Length(); ++i)
  {
    err = xhat->Update(alpha[i], *xbar, -alpha[i], *(sol[i]), 1.0);
    if (err != 0) { dserror("Update failed."); }
  }

  return xhat;
}

/*----------------------------------------------------------------------*/
/* Print summary of current iteration */
void NLNSOL::NlnOperatorNGmres::PrintIterSummary(const int iter,
    const double resnorm
    ) const
{
  // print only on one processor
  if (Comm().MyPID() == 0)
  {
    std::cout << std::setprecision(12)
              << "Finished NlnKrylov iteration " << iter
              << " with |f| = " << resnorm << std::endl;
  }
}
