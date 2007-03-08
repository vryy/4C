
#ifdef TRILINOS_PACKAGE

#include "fsi_nox_aitken.H"

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_GlobalData.H"

using namespace NOX;
using namespace NOX::LineSearch;

AitkenRelaxation::AitkenRelaxation(const Teuchos::RefCountPtr<NOX::Utils>& utils,
                                   Teuchos::ParameterList& params)
  : utils_(utils)
{
  Teuchos::ParameterList& p = params.sublist("Aitken");
  nu_ = p.get("Start nu", 0.0);
}

AitkenRelaxation::~AitkenRelaxation()
{
}

bool AitkenRelaxation::reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
                             Teuchos::ParameterList& params)
{
  Teuchos::ParameterList& p = params.sublist("Aitken");

  // do not reset the aitken factor
  //nu_ = p.get("Start nu", 0.0);

  if (!is_null(del_))
  {
    del_->init(1e20);
  }
  utils_ = gd->getUtils();
  return true;
}

bool AitkenRelaxation::compute(Abstract::Group& grp, double& step,
                               const Abstract::Vector& dir,
                               const Solver::Generic& s)
{
  if (utils_->isPrintType(NOX::Utils::InnerIteration))
  {
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n"
                  << "-- Aitken Line Search -- \n";
  }

  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  const NOX::Abstract::Vector& F = oldGrp.getF();

  if (is_null(del_))
  {
    del_  = F.clone(ShapeCopy);
    del2_ = F.clone(ShapeCopy);
    del_->init(1e20);
  }

  del2_->update(1,*del_,1,F);
  del_ ->update(-1,F);

  double top = del2_->innerProduct(*del_);
  double den = del2_->innerProduct(*del2_);

  nu_ = nu_ + (nu_ - 1.)*top/den;
  step = 1. - nu_;

  grp.computeX(oldGrp, dir, step);

#if 1
  // Why calculate F anew here? This is the second time in this
  // loop.
  // Do we need this in order to have an unrelaxed solution at the
  // time step end?
  grp.computeF();

  // is this reasonable at this point?
  double checkOrthogonality = fabs( grp.getF().innerProduct(dir) );

  if (utils_->isPrintType(Utils::InnerIteration)) {
    utils_->out() << setw(3) << "1" << ":";
    utils_->out() << " step = " << utils_->sciformat(step);
    utils_->out() << " orth = " << utils_->sciformat(checkOrthogonality);
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }
#else

  if (utils_->isPrintType(Utils::InnerIteration)) {
    utils_->out() << setw(3) << "1" << ":";
    utils_->out() << " step (omega) = " << utils_->sciformat(step);
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }

#endif

  return true;
}


#endif
