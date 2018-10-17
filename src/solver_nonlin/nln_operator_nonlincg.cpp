/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_nonlincg.cpp

\brief Nonlinear Conjugate Gradient

\level 3

\maintainer Matthias Mayr
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
      nlnprec_(Teuchos::null),
      betatype_(INPAR::NLNSOL::NONLINCG::beta_none),
      restartevery_(50)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::Setup()
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorNonlinCG::Setup");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  SetupLineSearch();
  SetupPreconditioner();

  // Determine the type of beta to be used
  const std::string betatype = MyGetParameter<std::string>("Nonlinear CG: Beta Type");
  if (betatype == "fletcherreeves")
    betatype_ = INPAR::NLNSOL::NONLINCG::beta_fletcherreeves;
  else if (betatype == "polakribiere")
    betatype_ = INPAR::NLNSOL::NONLINCG::beta_polakribiere;
  else if (betatype == "hestenesstiefel")
    betatype_ = INPAR::NLNSOL::NONLINCG::beta_hestenesstiefel;
  else
    dserror("'%s' is an unknown type for the parameter beta", betatype.c_str());

  restartevery_ = MyGetParameter<int>("Nonlinear CG: Restart Every Iterations");

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::SetupLineSearch()
{
  NLNSOL::LineSearchFactory linesearchfactory;
  linesearch_ =
      linesearchfactory.Create(Configuration(), MyGetParameter<std::string>("line search"));

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::SetupPreconditioner()
{
  const std::string opname = MyGetParameter<std::string>("Nonlinear CG: Nonlinear Preconditioner");

  NlnOperatorFactory nlnopfactory;
  nlnprec_ = nlnopfactory.Create(Configuration(), opname);
  nlnprec_->Init(Comm(), Configuration(), opname, NlnProblem(), BaciLinearSolver(), Nested() + 1);
  nlnprec_->Setup();

  return;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNonlinCG::ApplyInverse(
    const Epetra_MultiVector& f, Epetra_MultiVector& x) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnOperatorNonlinCG::ApplyInverse");
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

  // ---------------------------------------------------------------------------
  // some initializations
  // ---------------------------------------------------------------------------
  double alpha = 1.0;       // line search parameter
  double beta = 0.0;        // parameter for update of search direction
  double fnorm2 = 1.0e+12;  // L2-norm of residual
  int iter = 0;             // iteration counter
  bool suffdecr = false;    // flag for sufficient decrease of line search

  // ---------------------------------------------------------------------------
  // compute initial residual and apply the preconditioner once
  // ---------------------------------------------------------------------------
  // evaluate current residual
  Teuchos::RCP<Epetra_MultiVector> fnew = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  NlnProblem()->ComputeF(x, *fnew);

  // do an initial convergence check
  bool converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);
  PrintIterSummary(-1, fnorm2);

  // prepare vector for residual from previous iteration (needed for beta)
  Teuchos::RCP<Epetra_MultiVector> fold = Teuchos::rcp(new Epetra_MultiVector(*fnew));

  // Apply preconditioner once: s = M^{-1}*fnew
  Teuchos::RCP<Epetra_MultiVector> s = Teuchos::rcp(new Epetra_MultiVector(x));
  err = ApplyPreconditioner(*fnew, *s);
  if (err != 0)
  {
    dserror("ApplyPreconditioner() failed.");
  }
  s->Update(-1.0, x, 1.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  // prepare vector for preconditioned residual from previous iteration
  Teuchos::RCP<Epetra_MultiVector> sold = Teuchos::rcp(new Epetra_MultiVector(*s));

  // prepare vector for search direction and initialize to zero
  Teuchos::RCP<Epetra_MultiVector> p = Teuchos::rcp(new Epetra_MultiVector(s->Map(), 1, true));

  // do an initial convergence check
  converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);
  PrintIterSummary(iter, fnorm2);

  // ---------------------------------------------------------------------------
  // the nonlinear CG loop
  // ---------------------------------------------------------------------------
  while (ContinueIterations(iter, converged))
  {
    // Update search direction
    p->Update(1.0, *s, beta);

    // compute line search parameter alpha
    ComputeStepLength(x, *fnew, *p, fnorm2, alpha, suffdecr);

    // update solution
    err = x.Update(alpha, *p, 1.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }

    // update quantities from previous iteration
    fold->Update(1.0, *fnew, 0.0);
    sold->Update(1.0, *s, 0.0);

    // evaluate residual
    NlnProblem()->ComputeF(x, *fnew);
    converged = NlnProblem()->ConvergenceCheck(*fnew, fnorm2);

    // compute preconditioned search direction
    s->Update(1.0, x, 0.0);
    err = ApplyPreconditioner(*fnew, *s);
    if (err != 0)
    {
      dserror("ApplyPreconditioner() failed.");
    }
    s->Update(-1.0, x, 1.0);
    if (err != 0)
    {
      dserror("Update failed.");
    }

    ComputeBeta(beta, fnew, fold, s, sold);

    /* Account for possible restart:  We restart CG every #restartevery_
     * iterations or if beta <= 0.0. */
    if (((iter + 1) % restartevery_ == 0) or (beta <= 0.0))
    {
      beta = 0.0;

      if (getVerbLevel() > Teuchos::VERB_NONE)
      {
        *getOStream() << "   Restart " << Label() << " in iteration " << iter << std::endl;
      }
    }

    // finish current iteration
    ++iter;
    PrintIterSummary(iter, fnorm2);
  }

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

  // return error code
  return errorcode;
}

/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorNonlinCG::ApplyPreconditioner(
    const Epetra_MultiVector& f, Epetra_MultiVector& x) const
{
  if (nlnprec_.is_null())
    dserror(
        "Nonlinear preconditioner has not been initialized, yet. "
        "Has SetupPreconditioner() been called?");

  dsassert(f.Map().PointSameAs(x.Map()), "Maps do not match");

  return nlnprec_->ApplyInverse(f, x);
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::ComputeBeta(double& beta,
    Teuchos::RCP<const Epetra_MultiVector> fnew, Teuchos::RCP<const Epetra_MultiVector> fold,
    Teuchos::RCP<const Epetra_MultiVector> s, Teuchos::RCP<const Epetra_MultiVector> sold) const
{
  // compute beta based on user's choice
  switch (betatype_)
  {
    case INPAR::NLNSOL::NONLINCG::beta_none:
      dserror("Formula how to compute beta has not been set, yet.");
      break;
    case INPAR::NLNSOL::NONLINCG::beta_fletcherreeves:
      ComputeBetaFletcherReeves(beta, fnew, fold, s, sold);
      break;
    case INPAR::NLNSOL::NONLINCG::beta_polakribiere:
      ComputeBetaPolakRibiere(beta, fnew, fold, s, sold);
      break;
    case INPAR::NLNSOL::NONLINCG::beta_hestenesstiefel:
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
    Teuchos::RCP<const Epetra_MultiVector> fnew, Teuchos::RCP<const Epetra_MultiVector> fold,
    Teuchos::RCP<const Epetra_MultiVector> s, Teuchos::RCP<const Epetra_MultiVector> sold) const
{
  // compute numerator of beta
  fnew->Dot(*s, &beta);

  // compute denominator of beta
  double denominator = 1.0;
  fold->Dot(*sold, &denominator);

  // divide by the denominator
  beta /= denominator;

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::ComputeBetaPolakRibiere(double& beta,
    Teuchos::RCP<const Epetra_MultiVector> fnew, Teuchos::RCP<const Epetra_MultiVector> fold,
    Teuchos::RCP<const Epetra_MultiVector> s, Teuchos::RCP<const Epetra_MultiVector> sold) const
{
  // compute numerator of beta
  double num1 = 0.0;
  double num2 = 0.0;
  fnew->Dot(*s, &num1);
  fnew->Dot(*sold, &num2);

  // compute denominator of beta
  double denominator = 0.0;
  fold->Dot(*sold, &denominator);

  // compute beta
  beta = (num1 - num2) / denominator;

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::ComputeBetaHestenesStiefel(double& beta) const
{
  dserror("Not implemented, yet.");

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorNonlinCG::ComputeStepLength(const Epetra_MultiVector& x,
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
void NLNSOL::NlnOperatorNonlinCG::RebuildPrec()
{
  nlnprec_->RebuildPrec();

  return;
}
