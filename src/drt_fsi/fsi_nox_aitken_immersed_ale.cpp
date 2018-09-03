/*!----------------------------------------------------------------------
\file fsi_nox_aitken_immersed_ale.cpp

\brief class computing step length for AITKEN relaxation

\maintainer Andreas Rauch

\level 2

*----------------------------------------------------------------------*/

#include "fsi_nox_aitken_immersed_ale.H"
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

// AITKEN computed separate for immersed and ale part of interface
#include "../drt_immersed_problem/immersed_field_exchange_manager.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../linalg/linalg_mapextractor.H"

NOX::IMMERSED::AitkenRelaxationImmersedAle::AitkenRelaxationImmersedAle(
    const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params)
    : utils_(utils)
{
  Teuchos::ParameterList& p = params.sublist("Aitken");
  nu_ale_ = p.get("Start nu", 0.0);
  nu_immersed_ = nu_ale_;
  nu_ = nu_ale_;

  double maxstep = p.get("max step size", 0.0);
  if (maxstep > 0)
  {
    nu_ale_ = 1 - maxstep;
    nu_immersed_ = 1 - maxstep;
    nu_ = 1 - maxstep;
  }

  adapter_ = DRT::ImmersedFieldExchangeManager::Instance()->GetAdapter();
}


NOX::IMMERSED::AitkenRelaxationImmersedAle::~AitkenRelaxationImmersedAle() {}


bool NOX::IMMERSED::AitkenRelaxationImmersedAle::reset(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  Teuchos::ParameterList& p = params.sublist("Aitken");

  // We might want to constrain the step size of the first relaxation
  // in a new time step.
  double maxstep = p.get("max step size", 0.0);
  double minstep = p.get("min step size", 0.0);

  // check for maxstep
  if (maxstep > 0. && maxstep < 1. - nu_ale_) nu_ale_ = 1. - maxstep;

  if (maxstep > 0. && maxstep < 1. - nu_immersed_) nu_immersed_ = 1. - maxstep;

  if (maxstep > 0. && maxstep < 1. - nu_) nu_ = 1. - maxstep;

  // check for minstep
  if (minstep < 0. && minstep > 1. - nu_ale_) nu_ale_ = 1. - minstep;

  if (minstep < 0. && minstep > 1. - nu_immersed_) nu_immersed_ = 1. - minstep;

  if (minstep < 0. && minstep > 1. - nu_) nu_ = 1. - minstep;

  if (!is_null(del_))
  {
    del_->init(1e20);
  }
  utils_ = gd->getUtils();
  return true;
}


bool NOX::IMMERSED::AitkenRelaxationImmersedAle::compute(
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

  // turn off switch

  if (is_null(del_))
  {
    del_ = F.clone(ShapeCopy);
    del2_ = F.clone(ShapeCopy);
    del_->init(1.0e20);
    del2_->init(0.0);
  }

  // general update
  del2_->update(1, *del_, 1, F);
  del_->update(-1, F);

  double step_general = -1234.0;
  const double top = del2_->innerProduct(*del_);
  const double den = del2_->innerProduct(*del2_);

  nu_ = nu_ + (nu_ - 1.) * top / den;
  step_general = 1. - nu_;

  Teuchos::RCP<NOX::Epetra::Vector> epetra_del =
      Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector>(del_);
  if (epetra_del == Teuchos::null) dserror("cast failed");

  Teuchos::RCP<NOX::Epetra::Vector> epetra_del2 =
      Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector>(del2_);
  if (epetra_del == Teuchos::null) dserror("cast failed");

  // ale part
  double step_ale = -1234.0;
  Teuchos::RCP<Epetra_Vector> epetra_del2_ale =
      adapter_->CombinedInterface()->ExtractCondVector(epetra_del2->getEpetraVector());
  Teuchos::RCP<Epetra_Vector> epetra_del_ale =
      adapter_->CombinedInterface()->ExtractCondVector(epetra_del->getEpetraVector());

  double top_ale = -1234.0;
  epetra_del2_ale->Dot(*epetra_del_ale, &top_ale);
  double den_ale = -4321.0;
  epetra_del2_ale->Dot(*epetra_del2_ale, &den_ale);
  nu_ale_ = nu_ale_ + (nu_ale_ - 1.) * top_ale / den_ale;
  step_ale = 1. - nu_ale_;

  // immersed part
  double step_immersed = -1234.0;
  Teuchos::RCP<Epetra_Vector> epetra_del2_immersed =
      adapter_->CombinedInterface()->ExtractOtherVector(epetra_del2->getEpetraVector());
  Teuchos::RCP<Epetra_Vector> epetra_del_immersed =
      adapter_->CombinedInterface()->ExtractOtherVector(epetra_del->getEpetraVector());

  double top_immersed = -1234.0;
  epetra_del2_immersed->Dot(*epetra_del_immersed, &top_immersed);
  double den_immersed = -4321.0;
  epetra_del2_immersed->Dot(*epetra_del2_immersed, &den_immersed);

  nu_immersed_ = nu_immersed_ + (nu_immersed_ - 1.) * top_immersed / den_immersed;
  step_immersed = 1. - nu_immersed_;

  utils_->out() << "          RELAX ALE-FSI =      " << std::setprecision(13) << step_ale << "\n";
  utils_->out() << "          RELAX IMMERSED FSI = " << std::setprecision(13) << step_immersed
                << "\n";
  utils_->out() << "          OLD STEP =           " << std::setprecision(13) << step_general
                << "\n";
  utils_->out() << "          AVG. STEP =          " << std::setprecision(13)
                << ((step_ale + step_immersed) / 2.0) << "\n";

  // set step to 1.0 (only used as a dummy)
  step = (step_ale + step_immersed) / 2.0;

  // build new dir with different relaxation for immersed and ale part here and add this step with
  // dummy step size = 1.0 in computeX
  Teuchos::RCP<NOX::Epetra::Vector> combined_newdir =
      Teuchos::rcp(new NOX::Epetra::Vector(epetra_del->getEpetraVector()));
  combined_newdir->init(0.0);
  double norm = combined_newdir->norm();
  if (norm != 0.0) dserror("norm should be 0.0 here. Norm=%f", norm);

  // copy dir in plain dir
  Teuchos::RCP<NOX::Epetra::Vector> plain_dir =
      Teuchos::rcp(new NOX::Epetra::Vector(epetra_del->getEpetraVector()));
  plain_dir->init(0.0);
  norm = plain_dir->norm();
  if (norm != 0.0) dserror("norm should be 0.0 here. Norm=%f", norm);

  plain_dir->update(1.0, dir, 0.0);

  // extract and scale the ale and immersed parts of dir with their respective AITKEN factor
  Teuchos::RCP<Epetra_Vector> dir_ale =
      adapter_->CombinedInterface()->ExtractCondVector(plain_dir->getEpetraVector());
  dir_ale->Scale(step_ale);

  Teuchos::RCP<Epetra_Vector> dir_immersed =
      adapter_->CombinedInterface()->ExtractOtherVector(plain_dir->getEpetraVector());
  dir_immersed->Scale(step_immersed);

  // insert those two parts into combined direction vector
  adapter_->CombinedInterface()->InsertCondVector(
      dir_ale, Teuchos::rcpFromRef(combined_newdir->getEpetraVector()));
  adapter_->CombinedInterface()->InsertOtherVector(
      dir_immersed, Teuchos::rcpFromRef(combined_newdir->getEpetraVector()));

  // add combined new direction to x
  // grp.computeX(oldGrp, *combined_newdir, 1.0);
  // grp.computeX(oldGrp, dir, step_general);
  grp.computeX(oldGrp, dir, (step_ale + step_immersed) / 2.0);


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
    (*out) << count << " " << step_ale << " " << step_immersed << " " << fnorm << "\n";
    count += 1;
    out->flush();
  }

  return true;
}

// double NOX::IMMERSED::AitkenRelaxationImmersedAle::getOmega()
//{
//  return 1.-nu_immersed_;
//}
