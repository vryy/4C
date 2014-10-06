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
#include "nln_operator_base.H"
#include "nln_operator_factory.H"
#include "nln_operator_nonlincg.H"
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorNonlinCG::NlnOperatorNonlinCG()
: linesearch_(Teuchos::null),
  betatype_(NLNSOL::NlnOperatorNonlinCG::beta_none)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  SetupLineSearch();
  SetupPreconditioner();

  betatype_ = beta_fletcherreeves; // ToDo: Read this from input file/parameter list

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ =
      linesearchfactory.Create(Params().sublist("Nonlinear CG: Line Search"));

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::SetupPreconditioner()
{
  const Teuchos::ParameterList& precparams =
      Params().sublist("Nonlinear CG: Nonlinear Preconditioner");

  NlnOperatorFactory nlnopfactory;
  nlnprec_ = nlnopfactory.Create(precparams);
  nlnprec_->Init(Comm(), precparams, NlnProblem());
  nlnprec_->Setup();

  return;
}

/*----------------------------------------------------------------------------*/
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
  Teuchos::RCP<Epetra_MultiVector> fnew =
      Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->ComputeF(x, *fnew);

  // prepare vector for residual from previous iteration
  Teuchos::RCP<Epetra_MultiVector> fold =
      Teuchos::rcp(new Epetra_MultiVector(*fnew));

  // compute preconditioned search direction
  Teuchos::RCP<Epetra_MultiVector> p = Teuchos::rcp(new Epetra_MultiVector(x));
  ApplyPreconditioner(*fnew, x);
  p->Update(1.0, x, -1.0);

  bool converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);

  PrintIterSummary(iter, fnorm2);

  // ---------------------------------------------------------------------------
  // the nonlinear CG loop
  // ---------------------------------------------------------------------------
  while (ContinueIterations(iter, converged))
  {
    // compute line search parameter alpha
    alpha = ComputeStepLength(x, *p, fnorm2);

    // update solution
    err = x.Update(alpha, *p, 1.0);
    if (err != 0) { dserror("Update failed."); }

    // evaluate residual
    err = fold->Update(1.0, *fnew, 0.0);
    if (err != 0) { dserror("Update failed."); }
    NlnProblem()->ComputeF(x, *fnew);
    converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);

    // compute beta
    ComputeBeta(beta, *fnew, *fold);

    Teuchos::RCP<Epetra_MultiVector> pnew =
        Teuchos::rcp(new Epetra_MultiVector(p->Map(), true));

    // compute preconditioned search direction
    pnew->Update(1.0, x, 0.0);
    ApplyPreconditioner(*fnew, x);
    pnew->Update(1.0, x, -1.0);

    // Update
    err = p->Update(1.0, *pnew, beta);
    if (err != 0) { dserror("Update failed."); }

    ++iter;

    PrintIterSummary(iter, fnorm2);
  }

  CheckSuccessfulConvergence(iter, converged);

  // return error code
  return (not CheckSuccessfulConvergence(iter, converged));
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::ApplyPreconditioner(
    const Epetra_MultiVector& f, Epetra_MultiVector& x) const
{
  if (nlnprec_.is_null())
    dserror("Nonlinear preconditioner has not been initialized, yet. "
        "Has SetupPreconditioner() been called?");

  nlnprec_->ApplyInverse(f, x);

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::ComputeBeta(double& beta,
    const Epetra_MultiVector& fnew, const Epetra_MultiVector& fold) const
{
  // compute beta based on user's choice
  switch (betatype_)
  {
  case beta_none:
    dserror("Formula how to compute beta has not been set, yet.");
    break;
  case beta_fletcherreeves:
    ComputeBetaFletcherReeves(beta, fnew, fold);
    break;
  case beta_polakribiere:
    ComputeBetaPolakRibiere(beta);
    break;
  case beta_hestenesstiefel:
    ComputeBetaHestenesStiefel(beta);
    break;
  default:
    dserror("Unknown formula type to compute the parameter beta.");
    break;
  }

  // 'beta' has to be > 0.0, otherwise restart the procedure with 'beta = 0.0'
  beta = std::max(beta, 0.0);

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::ComputeBetaFletcherReeves(double& beta,
    const Epetra_MultiVector& fnew, const Epetra_MultiVector& fold) const
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
void NLNSOL::NlnOperatorNonlinCG::ComputeBetaPolakRibiere(double& beta) const
{
  dserror("Not implemented, yet.");

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::ComputeBetaHestenesStiefel(double& beta) const
{
  dserror("Not implemented, yet.");

  return;
}

/*----------------------------------------------------------------------------*/
const double NLNSOL::NlnOperatorNonlinCG::ComputeStepLength(
    const Epetra_MultiVector& x, const Epetra_MultiVector& inc,
    double fnorm2) const
{
  linesearch_->Init(NlnProblem(), Params().sublist("Nonlinear CG: Line Search"),
      x, inc, fnorm2);
  linesearch_->Setup();
  return linesearch_->ComputeLSParam();
}
