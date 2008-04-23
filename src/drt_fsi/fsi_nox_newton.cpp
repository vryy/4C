#ifdef CCADISCRET

#include "fsi_nox_newton.H"

#include <NOX_GlobalData.H>
#include <NOX_Utils.H>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::Newton::Newton(const Teuchos::RCP<NOX::GlobalData>& gd,
                         Teuchos::ParameterList& params)
  : NOX::Direction::Newton(gd,params),
    utils_(gd->getUtils()),
    desirednlnres_(0.),
    currentnlnres_(0.),
    plaintol_(1e-4),
    better_(0.1)
{
  // needed because NOX::Direction::Newton::Newton() does not call
  // NOX::FSI::Newton::reset()
  paramsPtr = &params;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool NOX::FSI::Newton::reset(const Teuchos::RCP<NOX::GlobalData>& gd,
                             Teuchos::ParameterList& params)
{
  paramsPtr = &params;
  return NOX::Direction::Newton::reset(gd,params);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool NOX::FSI::Newton::compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
                               const NOX::Solver::Generic& solver)
{
  //utils_->out() << "NOX::FSI::Newton::compute()\n";

  Teuchos::ParameterList& lsParams = paramsPtr->sublist("Newton").sublist("Linear Solver");
  lsParams.set<std::string>("Convergence Test","r0");

  // do adaptive linear solver tolerance

  // find largest current residual
  currentnlnres_ = 0;
  desirednlnres_ = 1;
  for (unsigned i=0; i<cresiduals_.size(); ++i)
  {
    double c = cresiduals_[i];
    double d = dresiduals_[i];
    if (currentnlnres_/desirednlnres_ < c/d)
    {
      currentnlnres_ = c;
      desirednlnres_ = d;
    }
  }

  utils_->out() << "                --- Aztec input   relative tolerance "
                << plaintol_
                << "\n";

  // heuristic tolerance calculation
  if (currentnlnres_*plaintol_ < desirednlnres_)
  {
    double tol = desirednlnres_*better_/currentnlnres_;
    if (tol>plaintol_)
    {
      lsParams.set<double>("Tolerance",tol);
      utils_->out() << "                *** Aztec adapted relative tolerance "
                    << tol
                    << "\n";
    }
  }
  else
  {
    lsParams.set<double>("Tolerance",plaintol_);
  }

  cresiduals_.clear();
  dresiduals_.clear();

  return NOX::Direction::Newton::compute(dir,grp,solver);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void NOX::FSI::Newton::Residual(double current, double desired)
{
  //utils_->out() << "NOX::FSI::Newton::Residual(" << current << "," << desired << ")\n";
  cresiduals_.push_back(current);
  dresiduals_.push_back(desired);
}


#endif
