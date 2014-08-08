/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_nonlincg.cpp

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

// baci
#include "linesearch_base.H"
#include "linesearch_factory.H"
#include "nln_operator_nonlincg.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperatorNonlinCG::NlnOperatorNonlinCG()
: linsolver_(Teuchos::null),
  linesearch_(Teuchos::null),
  maxiter_(1)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of the algorithm  / operator */
void NLNSOL::NlnOperatorNonlinCG::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  if (Params().isParameter("Nonlinear CG: Max Iter"))
    maxiter_ = Params().get<int>("Nonlinear CG: Max Iter");

  SetupLinearSolver();
  SetupLineSearch();

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of linear solver */
void NLNSOL::NlnOperatorNonlinCG::SetupLinearSolver()
{
  // get the solver number used for structural problems
  const int linsolvernumber = Params().get<int>("Nonlinear CG: Linear Solver");

  // check if the solver ID is valid
  if (linsolvernumber == (-1))
    dserror("No valid linear solver defined!");

  linsolver_ = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                               Comm(),
                                               DRT::Problem::Instance()->ErrorFile()->Handle()));

  return;
}

/*----------------------------------------------------------------------------*/
/* Setup of line search */
void NLNSOL::NlnOperatorNonlinCG::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ = linesearchfactory.Create(Params().sublist("Nonlinear CG: Line Search"));

  return;
}

/*----------------------------------------------------------------------------*/
/* Apply the preconditioner */
int NLNSOL::NlnOperatorNonlinCG::ApplyInverse(const Epetra_MultiVector& f,
    Epetra_MultiVector& x) const
{
  int err = 0;

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // ---------------------------------------------------------------------------
  // some initializations
  // ---------------------------------------------------------------------------
  double alpha = 1.0; // line search parameter
  double beta = 1.0; // parameter for update of search direction
  double fnorm2 = 1.0e+12;
  int iter = 0; // iteration counter

  // ---------------------------------------------------------------------------
  // compute initial search direction
  // ---------------------------------------------------------------------------
  // evaluate current residual
  Teuchos::RCP<Epetra_MultiVector> fnew = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->Evaluate(x, *fnew);

  // prepare vector for residual from previous iteration
  Teuchos::RCP<Epetra_MultiVector> fold = Teuchos::rcp(new Epetra_MultiVector(*fnew));

  // initial search direction
  Teuchos::RCP<Epetra_MultiVector> p = Teuchos::rcp(new Epetra_MultiVector(*fnew));

  int linsolve_error = linsolver_->Solve(NlnProblem()->GetJacobianOperator(), p, fold, true, false, Teuchos::null);
  if (linsolve_error != 0) { dserror("Linear solver failed."); }

  p->Scale(-1.0);

  bool converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);

  if (Params().get<bool>("Nonlinear CG: Print Iterations"))
    PrintIterSummary(iter, fnorm2);

  // ---------------------------------------------------------------------------
  // the nonlinear CG loop
  // ---------------------------------------------------------------------------
  while (iter < GetMaxIter() && not converged)
  {
    // compute line search parameter alpha
    linesearch_->Init(NlnProblem(), Params().sublist("Nonlinear CG: Line Search"), x, *p, fnorm2);
    linesearch_->Setup();
    alpha = linesearch_->ComputeLSParam();

    // update solution
    err = x.Update(alpha, *p, 1.0);
    if (err != 0) { dserror("Update failed."); }

    // evaluate residual
    err = fold->Update(1.0, *fnew, 0.0);
    if (err != 0) { dserror("Update failed."); }
    NlnProblem()->Evaluate(x, *fnew);
    converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);

    // compute beta
    ComputeBeta(beta, *fnew, *fold);

    Teuchos::RCP<Epetra_MultiVector> pnew = Teuchos::rcp(new Epetra_MultiVector(p->Map(), true));

    // do not refactor since matrix has not changed
    linsolve_error = linsolver_->Solve(NlnProblem()->GetJacobianOperator(), pnew, fnew, false, false, Teuchos::null);
    if (linsolve_error != 0) { dserror("Linear solver failed."); }

    // Update
    err = p->Update(-1.0, *pnew, beta);
    if (err != 0) { dserror("Update failed."); }

    ++iter;

    if (Params().get<bool>("Nonlinear CG: Print Iterations"))
      PrintIterSummary(iter, fnorm2);
  }

  // ---------------------------------------------------------------------------
  // check successful convergence
  // ---------------------------------------------------------------------------
  if (iter < GetMaxIter() && converged)
    return 0;
  else
  {
    if (IsSolver())
      dserror("Nonlinear CG did not converge in %d iterations.", iter);
  }

  // suppose non-convergence
  return 1;
}

/*----------------------------------------------------------------------------*/
/* Compute the parameter beta */
void NLNSOL::NlnOperatorNonlinCG::ComputeBeta(double& beta,
    const Epetra_MultiVector& fnew,
    const Epetra_MultiVector& fold
    ) const
{
  // ToDo We need a decision which beta is used
  ComputeBetaFletcherReeves(beta, fnew, fold);

//  ComputeBetaPolakRibiere(beta);
//  ComputeBetaHestenesStiefel(beta);

  // 'beta' has to be > 0.0, otherwise restart the procedure with 'beta = 0.0'
  beta = std::max(beta, 0.0);

  return;
}

/*----------------------------------------------------------------------------*/
/* Compute the parameter beta according to Fletcher-Reeves formula */
void NLNSOL::NlnOperatorNonlinCG::ComputeBetaFletcherReeves(double& beta,
  const Epetra_MultiVector& fnew,
  const Epetra_MultiVector& fold
  ) const
{
  // compute numerator of beta
  fnew.Dot(fnew, &beta);

  // compute denominator of beta
  double denominator = 1.0;
  fold.Dot(fold, &denominator);

  // divide by the denominator
  beta /= denominator;

  return;
}

/*----------------------------------------------------------------------------*/
/* Compute the parameter beta according to Polak-Ribiere formula */
void NLNSOL::NlnOperatorNonlinCG::ComputeBetaPolakRibiere(double& beta) const
{
  dserror("Not implemented, yet.");

  return;
}

/*----------------------------------------------------------------------------*/
/* Compute the parameter beta according to Hestenes-Stiefel formula */
void NLNSOL::NlnOperatorNonlinCG::ComputeBetaHestenesStiefel(double& beta) const
{
  dserror("Not implemented, yet.");

  return;
}
