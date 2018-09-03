/*----------------------------------------------------------------------------*/
/*!
\file linesearch_base.cpp

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
#include <Epetra_MultiVector.h>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "linesearch_base.H"
#include "nln_problem_base.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::LineSearchBase::LineSearchBase()
    : nlnproblem_(Teuchos::null),
      config_(Teuchos::null),
      listname_(""),
      xold_(Teuchos::null),
      fold_(Teuchos::null),
      inc_(Teuchos::null),
      resnormold_(0.0),
      suffdecrtype_(INPAR::NLNSOL::LINESEARCH::suffdecr_none),
      isinit_(false),
      issetup_(false)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBase::Init(Teuchos::RCP<const NLNSOL::NlnProblemBase> nlnproblem,
    Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config, const std::string listname,
    const Epetra_MultiVector& xold, const Epetra_MultiVector& fold, const Epetra_MultiVector& inc,
    const double resnormold)
{
  // We need to enforce calling Setup() after Init()
  issetup_ = false;

  // Initialize the member variables
  nlnproblem_ = nlnproblem;
  config_ = config;
  listname_ = listname;
  xold_ = Teuchos::rcp(&xold, false);
  fold_ = Teuchos::rcp(&fold, false);
  inc_ = Teuchos::rcp(&inc, false);
  resnormold_ = resnormold;

  suffdecrtype_ = INPAR::NLNSOL::LINESEARCH::suffdecr_armijo;
  //  suffdecrtype_ = INPAR::NLNSOL::LINESEARCH::suffdecr_aredpred;
  //  suffdecrtype_ = INPAR::NLNSOL::LINESEARCH::suffdecr_loose;

  // some sanity checks
  if (GetFNormOld() < 0.0)
    dserror(
        "Old residual norm 'resnormold_' = %f, but has to be greater than "
        "0.0!",
        resnormold_);

  const std::string verblevel = MyGetParameter<std::string>("line search: verbosity");
  setVerbLevel(NLNSOL::UTILS::TranslateVerbosityLevelToTeuchos(verblevel));

  if (getVerbLevel() > Teuchos::VERB_HIGH)
  {
    *getOStream() << "Parameter list passed to Line Search algorithm:" << std::endl;
    Configuration()->GetSubList(MyListName()).print(*getOStream());
  }

  // Init() has been called
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::LineSearchBase::IsSufficientDecrease(const double fnorm2, const double lsparam) const
{
  bool issufficientdecrease = false;

  switch (suffdecrtype_)
  {
    case INPAR::NLNSOL::LINESEARCH::suffdecr_armijo:
    {
      issufficientdecrease = SufficientDecreaseArmijo(fnorm2, lsparam);
      break;
    }
    case INPAR::NLNSOL::LINESEARCH::suffdecr_aredpred:
    {
      issufficientdecrease = SufficientDecreaseARedPRed(fnorm2);
      break;
    }
    case INPAR::NLNSOL::LINESEARCH::suffdecr_loose:
    {
      issufficientdecrease = SufficientDecreaseLoose(fnorm2);
      break;
    }
    case INPAR::NLNSOL::LINESEARCH::suffdecr_none:
    {
      dserror("Type of sufficient decrease condition hasn't been set.");
      break;
    }
    default:
    {
      dserror("Unknown type of sufficient decrease condition.");
      break;
    }
  }

  return issufficientdecrease;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::LineSearchBase::SufficientDecreaseArmijo(
    const double fnorm2, const double lsparam) const
{
  bool issufficientdecrease = false;

  const double alpha = 1.0e-4;  // as recommended in literature

  if (fnorm2 < (1.0 - alpha * lsparam) * resnormold_) issufficientdecrease = true;

  return issufficientdecrease;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::LineSearchBase::SufficientDecreaseARedPRed(const double fnorm2) const
{
  bool issufficientdecrease = false;

  const double alpha = 1.0e-4;  // as recommended in [Eisenstat (1996)]
  const double eta = 0.95;      // forcing term parameter

  if (fnorm2 < (1.0 - alpha * (1.0 - eta)) * resnormold_) issufficientdecrease = true;

  return issufficientdecrease;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::LineSearchBase::SufficientDecreaseLoose(const double fnorm2) const
{
  bool issufficientdecrease = false;

  const double kappa = 1.0;

  if (fnorm2 <= kappa * resnormold_) issufficientdecrease = true;

  return issufficientdecrease;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBase::Safeguard(double& lsparam, const double lsparamold) const
{
  lsparam = std::min(lsparam, 0.5 * lsparamold);
  lsparam = std::max(lsparam, 0.1 * lsparamold);

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBase::ComputeF(const Epetra_MultiVector& x, Epetra_MultiVector& f) const
{
  dsassert(x.Map().PointSameAs(f.Map()), "Maps do not match.");

  nlnproblem_->ComputeF(x, f);

  return;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::LineSearchBase::ConvergenceCheck(const Epetra_MultiVector& f, double& fnorm2) const
{
  return nlnproblem_->ConvergenceCheck(f, fnorm2);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> NLNSOL::LineSearchBase::Configuration() const
{
  // check if parameter list has already been set
  if (config_.is_null()) dserror("Configuration 'config_' has not been initialized, yet.");

  return config_;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBase::ComputeLSParam(double& lsparam) const
{
  bool suffdecr = false;
  ComputeLSParam(lsparam, suffdecr);

  return;
}
