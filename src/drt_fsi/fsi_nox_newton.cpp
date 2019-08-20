/*----------------------------------------------------------------------*/
/*! \file

\brief NOX Newton direction with adaptive linear solver tolerance for FSI

\level 2

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#include "fsi_nox_newton.H"

#include <NOX_GlobalData.H>
#include <NOX_Utils.H>

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::Newton::Newton(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
    : NOX::Direction::Newton(gd, params),
      utils_(gd->getUtils()),
      desirednlnres_(0.),
      currentnlnres_(0.),
      plaintol_(1e-4),
      better_(0.1),
      verbosity_(INPAR::FSI::verbosity_full)
{
  // needed because NOX::Direction::Newton::Newton() does not call
  // NOX::FSI::Newton::reset()
  paramsPtr = &params;

  // adaptive tolerance settings for linear solver
  Teuchos::ParameterList& lsParams = paramsPtr->sublist("Newton").sublist("Linear Solver");
  plaintol_ = lsParams.get<double>("base tolerance");             // relative tolerance
  better_ = lsParams.get<double>("adaptive distance");            // adaptive distance
  verbosity_ = lsParams.get<INPAR::FSI::Verbosity>("verbosity");  // verbosity level
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool NOX::FSI::Newton::reset(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  paramsPtr = &params;

  Teuchos::ParameterList& lsParams = paramsPtr->sublist("Newton").sublist("Linear Solver");
  plaintol_ = lsParams.get<double>("base tolerance");

  return NOX::Direction::Newton::reset(gd, params);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool NOX::FSI::Newton::compute(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, const NOX::Solver::Generic& solver)
{
  TEUCHOS_FUNC_TIME_MONITOR("NOX::FSI::Newton::compute");

  Teuchos::ParameterList& lsParams = paramsPtr->sublist("Newton").sublist("Linear Solver");
  lsParams.set<std::string>("Convergence Test", "r0");

  // do adaptive linear solver tolerance

  // find largest current residual
  currentnlnres_ = 0;
  desirednlnres_ = 1;
  if (better_ > 0)
  {
    for (unsigned i = 0; i < cresiduals_.size(); ++i)
    {
      double c = cresiduals_[i];
      double d = dresiduals_[i];
      if (currentnlnres_ / desirednlnres_ < c / d)
      {
        currentnlnres_ = c;
        desirednlnres_ = d;
      }
    }
  }

  if (verbosity_ >= INPAR::FSI::verbosity_medium)
  {
    utils_->out() << "                --- Aztec input   relative tolerance " << plaintol_ << "\n";
  }

  // Always reset. To be sure.
  lsParams.set<double>("Tolerance", plaintol_);

  // heuristic tolerance calculation
  if (better_ > 0)
  {
    if (currentnlnres_ * plaintol_ < desirednlnres_)
    {
      double tol = desirednlnres_ * better_ / currentnlnres_;
      if (tol > plaintol_)
      {
        lsParams.set<double>("Tolerance", tol);
        if (verbosity_ >= INPAR::FSI::verbosity_medium)
        {
          utils_->out() << "                *** Aztec adapted relative tolerance " << tol << "\n";
        }
      }
    }
  }

  cresiduals_.clear();
  dresiduals_.clear();

  return NOX::Direction::Newton::compute(dir, grp, solver);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void NOX::FSI::Newton::Residual(double current, double desired)
{
  // utils_->out() << "NOX::FSI::Newton::Residual(" << current << "," << desired << ")\n";
  cresiduals_.push_back(current);
  dresiduals_.push_back(desired);
}
