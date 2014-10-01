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

// standard
#include <iostream>

// Epetra
#include <Epetra_MultiVector.h>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "linesearch_base.H"
#include "nln_problem.H"

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::LineSearchBase::LineSearchBase()
 : nlnproblem_(Teuchos::null),
   params_(Teuchos::null),
   xold_(Teuchos::null),
   inc_(Teuchos::null),
   resnormold_(0.0),
   isinit_(false),
   issetup_(false)
{
  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBase::Init(
    Teuchos::RCP<const NLNSOL::NlnProblem> nlnproblem,
    const Teuchos::ParameterList& params,
    const Epetra_MultiVector& xold,
    const Epetra_MultiVector& inc,
    const double resnormold)
{
  // We need to enforce calling Setup() after Init()
  issetup_ = false;

  // Initialize the member variables
  nlnproblem_ = nlnproblem;
  params_ = Teuchos::rcp(&params, false);
  xold_ = Teuchos::rcp(&xold, false);
  inc_ = Teuchos::rcp(&inc, false);
  resnormold_ = resnormold;

  // some sanity checks
  if (GetFNormOld() < 0.0)
    dserror("Old residual norm 'resnormold_' = %f, but has to be greater than "
        "0.0!", resnormold_);

  // Init() has been called
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::LineSearchBase::IsSufficientDecrease(const double fnorm2,
    const double lsparam
    ) const
{
  bool issufficientdecrease = false;

  const double alpha = 1.0e-4; // as recommended in literature

  if (fnorm2 < (1.0 - alpha * lsparam) * resnormold_)
    issufficientdecrease = true;

  return issufficientdecrease;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBase::Safeguard(double& lsparam,
    const double lsparamold
    ) const
{
  lsparam = std::min(lsparam, 0.5*lsparamold);
  lsparam = std::max(lsparam, 0.1*lsparamold);

  return;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::LineSearchBase::ComputeF(const Epetra_MultiVector& x,
    Epetra_MultiVector& f
    ) const
{
  dsassert(x.Map().PointSameAs(f.Map()), "Maps do not match.");

  nlnproblem_->ComputeF(x, f);

  return;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::LineSearchBase::ConvergenceCheck(const Epetra_MultiVector& f,
    double& fnorm2
    ) const
{
  return nlnproblem_->ConvergenceCheck(f, fnorm2);
}

/*----------------------------------------------------------------------------*/
const Teuchos::ParameterList& NLNSOL::LineSearchBase::Params() const
{
  // check if parameter list has already been set
  if (params_.is_null())
    dserror("Parameter list 'params_' has not been initialized, yet.");

  return *params_;
}
