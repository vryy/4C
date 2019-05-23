/*!----------------------------------------------------------------------

\brief class computing step length for AITKEN relaxation

\maintainer Matthias Mayr

\level 1

*----------------------------------------------------------------------*/

#include "fsi_nox_aitken.H"
#include "fsi_utils.H"

#include <NOX_Common.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Solver_Generic.H>
#include <Teuchos_ParameterList.hpp>
#include <NOX_GlobalData.H>

#include <Epetra_Vector.h>
#include <Epetra_Comm.h>
#include <NOX_Epetra_Vector.H>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

NOX::FSI::AitkenRelaxation::AitkenRelaxation(
    const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params)
    : utils_(utils)
{
  Teuchos::ParameterList& p = params.sublist("Aitken");
  nu_ = p.get("Start nu", 0.0);

  maxstep_ = p.get("max step size", 0.0);
  if (maxstep_ > 0.0) nu_ = 1.0 - maxstep_;

  minstep_ = p.get("min step size", -1.0);

  restart_ = p.get<int>("restart", 0);

  restart_omega_ = p.get<double>("restart_omega", 0.0);
}


NOX::FSI::AitkenRelaxation::~AitkenRelaxation() {}


bool NOX::FSI::AitkenRelaxation::reset(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  // We might want to constrain the step size of the first relaxation
  // in a new time step.
  if (maxstep_ > 0.0 && maxstep_ < 1.0 - nu_) nu_ = 1. - maxstep_;

  if (minstep_ > 1.0 - nu_) nu_ = 1.0 - minstep_;

  if (!is_null(del_))
  {
    del_->init(1e20);
  }
  utils_ = gd->getUtils();
  return true;
}


bool NOX::FSI::AitkenRelaxation::compute(
    Abstract::Group& grp, double& step, const Abstract::Vector& dir, const Solver::Generic& s)
{
  if (utils_->isPrintType(NOX::Utils::InnerIteration))
  {
    utils_->out() << "\n"
                  << NOX::Utils::fill(72) << "\n"
                  << "-- Aitken Line Search -- \n";
  }

  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  const NOX::Abstract::Vector& F = oldGrp.getF();


  if (is_null(del_))
  {
    del_ = F.clone(ShapeCopy);
    del2_ = F.clone(ShapeCopy);
    del_->init(1.0e20);
    del2_->init(0.0);
  }

  del2_->update(1, *del_, 1, F);
  del_->update(-1, F);

  const double top = del2_->innerProduct(*del_);
  const double den = del2_->innerProduct(*del2_);

  // in case of a restart, we use the omega that
  // was written as restart output
  if (not restart_)
    nu_ = nu_ + (nu_ - 1.) * top / den;
  else
    nu_ = 1.0 - restart_omega_;

  // check constraints for step size
  if (maxstep_ > 0.0 && maxstep_ < 1.0 - nu_) nu_ = 1. - maxstep_;

  if (minstep_ > 1.0 - nu_) nu_ = 1.0 - minstep_;

  // calc step
  step = 1. - nu_;

  utils_->out() << "          RELAX = " << std::setw(5) << step << "\n";

  grp.computeX(oldGrp, dir, step);

  // Calculate F anew here. This results in another FSI loop. However
  // the group will store the result, so it will be reused until the
  // group's x is changed again. We do not waste anything.
  grp.computeF();

  // is this reasonable at this point?
  double checkOrthogonality = fabs(grp.getF().innerProduct(dir));

  if (utils_->isPrintType(Utils::InnerIteration))
  {
    utils_->out() << std::setw(3) << "1"
                  << ":";
    utils_->out() << " step = " << utils_->sciformat(step);
    utils_->out() << " orth = " << utils_->sciformat(checkOrthogonality);
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  // write omega
  double fnorm = grp.getF().norm();
  if (dynamic_cast<const NOX::Epetra::Vector&>(F).getEpetraVector().Comm().MyPID() == 0)
  {
    static int count;
    static std::ofstream* out;
    if (out == NULL)
    {
      std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
      s.append(".omega");
      out = new std::ofstream(s.c_str());
    }
    (*out) << count << " " << step << " " << fnorm << "\n";
    count += 1;
    out->flush();
  }

  // reset restart flag
  restart_ = false;

  return true;
}


double NOX::FSI::AitkenRelaxation::GetOmega() { return 1. - nu_; }
