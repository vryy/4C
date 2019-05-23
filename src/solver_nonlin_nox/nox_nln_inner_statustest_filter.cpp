/*----------------------------------------------------------------------------*/
/*!

\brief Inner status test class for constraint problems. Filter
       techniques are based on ideas from multi-objective optimization:

       - Control of the two distinct goals of minimization of the objective
         function and satisfaction of the constraints.

       - Unlike merit functions, filter methods keep these two goals separate

\maintainer Anh-Tu Vuong


\level 3

*/
/*----------------------------------------------------------------------------*/

#include "nox_nln_inner_statustest_filter.H"
#include "nox_nln_inner_statustest_interface_required.H"
#include "nox_nln_meritfunction_lagrangian.H"
#include "nox_nln_linesearch_generic.H"
#include "nox_nln_solver_linesearchbased.H"
#include "nox_nln_statustest_normf.H"
#include "nox_nln_statustest_activeset.H"
#include "nox_nln_group.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_Time.hpp>
#include <NOX_MeritFunction_Generic.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Direction_Generic.H>
#include <NOX_Epetra_Vector.H>
#include <Epetra_Vector.h>

#include <fenv.h>

/*----------------------------------------------------------------------------*/
// definition and initialization of static members
unsigned NOX::NLN::INNER::StatusTest::Filter::Point::num_coords_(0);
unsigned NOX::NLN::INNER::StatusTest::Filter::Point::num_obj_coords_(0);
double NOX::NLN::INNER::StatusTest::Filter::Point::gamma_obj_(0.0);
double NOX::NLN::INNER::StatusTest::Filter::Point::gamma_theta_(0.0);
std::vector<bool> NOX::NLN::INNER::StatusTest::Filter::Point::isvalid_scaling_;
LINALG::SerialDenseVector NOX::NLN::INNER::StatusTest::Filter::Point::scale_;
LINALG::SerialDenseVector NOX::NLN::INNER::StatusTest::Filter::Point::weights_;
LINALG::SerialDenseVector NOX::NLN::INNER::StatusTest::Filter::Point::global_scaled_max_thetas_;
double NOX::NLN::INNER::StatusTest::Filter::Point::global_init_max_theta_scale_;
std::set<Teuchos::RCP<NOX::NLN::INNER::StatusTest::Filter::Point>,
    NOX::NLN::INNER::StatusTest::rcp_comp<NOX::NLN::INNER::StatusTest::Filter::Point>>
    NOX::NLN::INNER::StatusTest::Filter::Point::filter_point_register_;

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::Filter::Filter(const FilterParams& fparams, const NOX::Utils& utils)
    : status_(status_unevaluated),
      theta_(fparams.infeasibility_vec_),
      curr_points_(plain_const_point_pair(Teuchos::null, Teuchos::null)),
      curr_fpoints_(plain_point_pair(Teuchos::null, Teuchos::null)),
      filter_(),
      non_dominated_filter_points_(),
      backup_(),
      soc_(SOCBase::create(*this, fparams.use_soc_, fparams.soc_type_)),
      blocking_ptr_(Blocking::create(*this, fparams)),
      blocking_(*blocking_ptr_),
      gamma_alpha_(fparams.gamma_alpha_),
      amin_obj_(1.0),
      amin_theta_(1.0),
      amin_ftype_(1.0),
      amin_(1.0),
      sf_(fparams.sf_),
      st_(fparams.st_),
      model_lin_terms_(),
      model_mixed_terms_(),
      armijo_test_(fparams.armijo_),
      is_ftype_step_(false),
      theta_min_ftype_(1.0),
      filter_status_(FilterStatusType::unevaluated),
      utils_(utils)
{
  Point::resetStaticMembers(1, theta_.number_, fparams.weight_objective_func_,
      fparams.weight_infeasibility_func_, blocking_.init_max_theta_scaling_);

  Point::setMarginSafetyFactors();

  model_lin_terms_.Resize(Point::num_coords_);
  model_mixed_terms_.Resize(Point::num_coords_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Reset()
{
  filter_.clear();
  blocking_.filter_iterates_.clear();
  Point::resetStaticMembers();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::setMarginSafetyFactors()
{
  gamma_obj_ = std::min<double>(
      1.0e-6, 1.0 / (2.0 * std::sqrt<double>(static_cast<double>(num_coords_)).real()));

  gamma_theta_ = gamma_obj_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::clearFilterPointRegister()
{
  filter_point_register_.clear();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::resetStaticMembers()
{
  if (isvalid_scaling_.size() < num_coords_)
    isvalid_scaling_.resize(num_coords_, false);
  else
    std::fill(isvalid_scaling_.begin(), isvalid_scaling_.end(), false);

  scale_.LightResize(num_coords_);

  // initialize the scale_ vector
  std::fill(scale_.A(), scale_.A() + num_coords_, 1.0);

  const unsigned num_theta = num_coords_ - num_obj_coords_;
  global_scaled_max_thetas_.LightResize(num_theta);

  std::fill(global_scaled_max_thetas_.A(), global_scaled_max_thetas_.A() + num_theta, 0.0);

  clearFilterPointRegister();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::resetStaticMembers(const unsigned num_obj_coords,
    const unsigned num_theta_coords, const double weight_objective_func,
    const double weight_infeasibility_func, const double init_max_theta_scale)
{
  num_obj_coords_ = num_obj_coords;
  num_coords_ = num_theta_coords + num_obj_coords;

  weights_.LightResize(num_coords_);

  // initialize the weights vector
  std::fill(weights_.A(), weights_.A() + num_obj_coords_, weight_objective_func);
  std::fill(weights_.A() + num_obj_coords_, weights_.A() + num_coords_, weight_infeasibility_func);

  global_init_max_theta_scale_ = init_max_theta_scale;

  resetStaticMembers();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::reinitFilter(
    plain_point_set& filter, const Infeasibility& infeasibility_func, const double& downscale_fac)
{
  filter.clear();
  clearFilterPointRegister();

  scaleMaxThetaValues(downscale_fac);

  Teuchos::RCP<Point> point_ptr = Teuchos::rcp(new Point);
  Point& point = *point_ptr;

  const unsigned num_thetas = global_scaled_max_thetas_.Length();
  std::copy(global_scaled_max_thetas_.A(), global_scaled_max_thetas_.A() + num_thetas,
      point.coords_.A() + num_obj_coords_);

  point.max_theta_id_ = infeasibility_func.findMaxThetaId(point.A() + Point::num_obj_coords_);

  Teuchos::RCP<Point> fpoint = makeFilterPoint(point, false);
  filter.push_back(fpoint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Filter::Point>
NOX::NLN::INNER::StatusTest::Filter::Point::create(const NOX::MeritFunction::Generic& merit_func,
    const Infeasibility& infeasibility_func, const NOX::Abstract::Group& grp)
{
  Teuchos::RCP<Point> point_ptr = Teuchos::rcp(new Point);
  Point& point = *point_ptr;

  point(0) = merit_func.computef(grp);
  infeasibility_func.computef(point.A() + Point::num_obj_coords_, grp);

  point.max_theta_id_ = infeasibility_func.findMaxThetaId(point.A() + Point::num_obj_coords_);

  return point_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Filter::Point>
NOX::NLN::INNER::StatusTest::Filter::Point::makeFilterPoint(const Point& p, const bool do_scaling)
{
  Teuchos::RCP<Point> fp_ptr = Teuchos::rcp(new Point(p));
  Point& fp = *fp_ptr;

  if (not fp.is_filter_point_)
  {
    if (do_scaling) fp.scale();
    fp.setNorm();
    fp.setMargin();
    fp.is_filter_point_ = true;
  }

  addFilterPointToRegister(fp_ptr);
  return fp_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::addFilterPointToRegister(
    const Teuchos::RCP<Point>& fp_ptr)
{
  for (const Teuchos::RCP<Point>& reg_fp_ptr : filter_point_register_)
  {
    // erase all old filter points which are no longer in use
    if (reg_fp_ptr.strong_count() == 1) filter_point_register_.erase(reg_fp_ptr);
  }

  filter_point_register_.insert(fp_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::scaleCoordinateOfAllRegisteredFilterPoints(
    const int id)
{
  if (isvalid_scaling_[id]) return;

  for (const Teuchos::RCP<Point>& reg_fp_ptr : filter_point_register_)
  {
    reg_fp_ptr->coords_(id) *= scale_(id);
    reg_fp_ptr->setNorm();
    reg_fp_ptr->setMargin();
  }

  isvalid_scaling_[id] = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::scaleMaxThetaValues(const double& fac)
{
  const unsigned length = static_cast<unsigned>(global_scaled_max_thetas_.Length());
  for (unsigned i = 0; i < length; ++i) global_scaled_max_thetas_(i) *= fac;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::setInitialScaledMaxThetaValue(
    const int id, const double& val)
{
  if (isvalid_scaling_[id]) return;

  const int theta_id = id - static_cast<int>(num_obj_coords_);
  if (theta_id < 0) return;

  // we take the maximal theta values and multiply it by two to specify an upper
  // bound
  global_scaled_max_thetas_(theta_id) = global_init_max_theta_scale_ * val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::scale()
{
  for (int i = 0; i < coords_.Length(); ++i)
  {
    if (coords_(i) != 0.0)
    {
      if (not isvalid_scaling_[i])
      {
        scale_(i) = weights_(i) / std::abs(coords_(i));
        setInitialScaledMaxThetaValue(i, coords_(i) * scale_(i));
        scaleCoordinateOfAllRegisteredFilterPoints(i);
      }

      coords_(i) *= scale_(i);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::setMargin()
{
  const double max_theta = maxTheta();

  for (unsigned i = 0; i < num_obj_coords_; ++i)
  {
    margin_(i) = gamma_obj_ * max_theta;
  }

  for (unsigned i = num_obj_coords_; i < num_coords_; ++i)
  {
    margin_(i) = gamma_theta_ * max_theta;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::INNER::StatusTest::Filter::Point::IsFeasible(const double tol) const
{
  if (is_feasible_) return true;

  if (not is_filter_point_) dserror("The point is expected to be a filter point.");

  const double used_tol = std::max(0.0, tol);
  for (unsigned i = num_obj_coords_; i < num_coords_; ++i)
    if (coords_(i) <= used_tol * scale_(i))
    {
      is_feasible_ = true;
      return true;
    }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::INNER::StatusTest::Filter::Point::IsSufficientlyReducedComparedToMaxTheta(
    const double& red_fac) const
{
  if (not is_filter_point_) dserror("This routine is only meant for filter points!");

  unsigned j = 0;
  for (unsigned i = num_obj_coords_; i < num_coords_; ++i, ++j)
    if (coords_(i) > red_fac * global_scaled_max_thetas_(j)) return false;

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::setNorm() { norm_ = coords_.Norm2(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::INNER::StatusTest::Filter::Point::print(
    std::ostream& stream, int par_indent_length, const NOX::Utils* u) const
{
  std::string par_indent;
  par_indent.assign(par_indent_length, ' ');

  stream << par_indent;
  if (is_filter_point_)
  {
    stream << "Filter-Point -- { ";
  }
  else
  {
    stream << "Point -- { ";
  }

  for (unsigned i = 0; i < num_coords_; ++i)
  {
    if (i != 0) stream << ", ";
    stream << NOX::Utils::sciformat(coords_(i), OUTPUT_PRECISION);
  }
  stream << " } with norm = " << NOX::Utils::sciformat(norm_, OUTPUT_PRECISION) << ";\n";

  if (!u or u->isPrintType(NOX::Utils::Details))
  {
    stream << par_indent << "margin = { ";
    for (unsigned i = 0; i < num_coords_; ++i)
    {
      if (i != 0) stream << ", ";
      stream << NOX::Utils::sciformat(margin_(i), OUTPUT_PRECISION);
    }
    stream << " };\n";

    stream << par_indent << "MaxTheta = { "
           << "id = " << max_theta_id_ << ", value = " << maxTheta()
           << ", scale = " << NOX::Utils::sciformat(scaleOfMaxTheta(), OUTPUT_PRECISION) << " };\n";

    stream << par_indent << "IsFeasible = " << (is_feasible_ ? "TRUE" : "FALSE") << "\n";
  }


  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Infeasibility::computef(
    double* theta_values, const NOX::Abstract::Group& grp) const
{
  plain_merit_func_set::const_iterator cit = vector_.begin();

  for (unsigned i = 0; cit != vector_.end(); ++cit, ++i)
  {
    const NOX::MeritFunction::Generic& func = **cit;
    theta_values[i] = func.computef(grp);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned NOX::NLN::INNER::StatusTest::Filter::Infeasibility::findMaxThetaId(
    double* theta_values) const
{
  unsigned max_theta_id = 0;
  unsigned i = 1;

  while (i < number_)
  {
    if (theta_values[max_theta_id] < theta_values[i]) max_theta_id = i;

    ++i;
  }

  return max_theta_id;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::InitPoints(const Interface::Required& interface,
    const NOX::Solver::Generic& solver, const NOX::Abstract::Group& grp)
{
  const int iter_newton = solver.getNumIterations();
  const NOX::MeritFunction::Generic& merit_func = interface.GetMeritFunction();

  switch (iter_newton)
  {
    case 0:
    {
      dserror("Not supposed to be called for Newton iteration zero!");
      break;
    }
    // --------------------------------------------------------------------
    // compute the point coordinates of the first reference state
    // which is accepted by default and clear the old filter
    // --------------------------------------------------------------------
    case 1:
    {
      Reset();
      curr_points_.first = Point::create(merit_func, theta_, grp);
      curr_fpoints_.first = Point::makeFilterPoint(*curr_points_.first, false);

      // Note that the value is currently set to 1.0e-4, if the predictor
      // step does not cause any penetration.
      theta_min_ftype_ = 1.0e-4 * std::max(1.0, curr_points_.first->maxTheta());

      break;
    }
    // --------------------------------------------------------------------
    // Move the previously accepted point to the first position at the very
    // beginning of each Newton step (except for the very first Newton step)
    // --------------------------------------------------------------------
    default:
    {
      // set accepted trial point at the first position
      if (not curr_points_.second.is_null()) curr_points_.first = curr_points_.second;
      if (not curr_fpoints_.second.is_null()) curr_fpoints_.first = curr_fpoints_.second;

      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::Filter::CheckStatus(
    const Interface::Required& interface, const NOX::Solver::Generic& solver,
    const NOX::Abstract::Group& grp, NOX::StatusTest::CheckType checkType)
{
  const NOX::NLN::LineSearch::Generic* linesearch_ptr =
      dynamic_cast<const NOX::NLN::LineSearch::Generic*>(&interface);
  if (not linesearch_ptr) dserror("Dynamic cast failed!");
  const NOX::NLN::LineSearch::Generic& linesearch = *linesearch_ptr;

  // do stuff at the beginning of a line search call
  const int iter_ls = interface.GetNumIterations();
  if (iter_ls == 0)
  {
    const NOX::Abstract::Vector& dir = linesearch.GetSearchDirection();

    InitPoints(interface, solver, grp);

    // setup the linear and quadratic model terms in the beginning
    SetupModelTerms(dir, grp, interface);

    // compute the minimal step length estimates
    ComputeMinimalStepLengthEstimates();

    // setup armijo test
    armijo_test_->CheckStatus(interface, solver, grp, checkType);

    // create a back-up of the last accepted state
    backup_.create(grp, dir);

    status_ = status_unevaluated;

    return status_;
  }
  else if (checkType == NOX::StatusTest::None)
  {
    status_ = status_unevaluated;
    filter_status_ = FilterStatusType::unevaluated;

    return status_;
  }

  ExecuteCheckStatus(linesearch, solver, grp, checkType);

  return PostCheckStatus(linesearch, solver, grp, checkType);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::SetTrialPoint(
    const NOX::MeritFunction::Generic& merit_func, const NOX::Abstract::Group& grp)
{
  curr_points_.second = Point::create(merit_func, theta_, grp);
  curr_fpoints_.second = Point::makeFilterPoint(*curr_points_.second, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::ExecuteCheckStatus(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Solver::Generic& solver,
    const NOX::Abstract::Group& grp, NOX::StatusTest::CheckType checkType)
{
  const NOX::MeritFunction::Generic& merit_func = linesearch.GetMeritFunction();

  // reset the f-type flag
  is_ftype_step_ = false;

  // compute the new point coordinates of the trial point
  SetTrialPoint(merit_func, grp);

  // trial filter point
  const Point& trial_fp = *curr_fpoints_.second;
  filter_status_ = AcceptabilityCheck(trial_fp);

  // get current step length
  const double step = linesearch.GetStepLength();

  switch (filter_status_)
  {
    // if the current trial point is not in the taboo region, we will check a
    // 2-nd criterion
    case FilterStatusType::passed_point_by_point:
    {
      // ------------------------------------------
      // Final F-Type check
      // ------------------------------------------
      if (CheckFTypeSwitchingCondition(step))
      {
        is_ftype_step_ = true;
        status_ = armijo_test_->CheckStatus(linesearch, solver, grp, checkType);
      }
      // ------------------------------------------
      // Final filter check
      // ------------------------------------------
      else
      {
        status_ = SufficientReductionCheck(trial_fp);


        if (status_ == status_converged and trial_fp.norm_ >= 0.0 and
            not trial_fp.IsFeasible(GetConstraintTolerance(solver)))
          AugmentFilter();
      }

      if (status_ != status_converged) status_ = IsAdmissibleStep(solver, step);

      break;
    }
    case FilterStatusType::rejected:
    {
      status_ = IsAdmissibleStep(solver, step);
      blocking_.Check(linesearch, solver, grp, trial_fp);

      break;
    }
    default:
    {
      dserror(
          "Unexpected filter status type: "
          "The filer %s.",
          FilterStatus2String(filter_status_).c_str());
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::GetConstraintTolerance(
    const NOX::Solver::Generic& solver) const
{
  const NOX::NLN::Solver::LineSearchBased* ls_solver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&solver);
  if (not ls_solver) dserror("The given non-linear solver is not line search based!");

  NOX::StatusTest::Generic* normf_test =
      ls_solver->GetOuterStatusTestWithQuantity<NOX::NLN::StatusTest::NormF>(
          NOX::NLN::StatusTest::quantity_contact_normal);

  if (normf_test)
  {
    const double true_tol = dynamic_cast<NOX::NLN::StatusTest::NormF&>(*normf_test)
                                .GetTrueTolerance(NOX::NLN::StatusTest::quantity_contact_normal);
    if (true_tol < 0.0) dserror("Something went wrong!");

    utils_.out(NOX::Utils::Debug) << __LINE__ << " -- " << __PRETTY_FUNCTION__ << "\n"
                                  << "True tolerance of the CONTACT-NORMAL-F-NORM test is "
                                  << true_tol << std::endl;

    return true_tol;
  }

  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Blocking::Check(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Solver::Generic& solver,
    const NOX::Abstract::Group& grp, const Point& rejected_fp)
{
  // sanity check
  if (filter_.filter_status_ != FilterStatusType::rejected) return;

  // Would the filter point be theoretically accepted by the inner filter test
  StatusType trial_status = status_unevaluated;
  if (filter_.CheckFTypeSwitchingCondition(linesearch.GetStepLength()))
    trial_status =
        filter_.armijo_test_->CheckStatus(linesearch, solver, grp, NOX::StatusTest::Complete);
  else
    trial_status = filter_.SufficientReductionCheck(rejected_fp);

  /* If the filter rejected the trial point but the inner test would accept it
   * we have an indicator for a blocking filter. This can happen due to old
   * historic information which is not reliable for the current neighborhood. */
  if (trial_status == status_converged)
  {
    const int newton_iter = solver.getNumIterations();
    AddFilterIterate(newton_iter, rejected_fp);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Blocking::AddFilterIterate(
    const int newton_iter, const Point& rejected_fp)
{
  if (newton_iter <= 0) return;

  const bool is_sufficiently_red =
      rejected_fp.IsSufficientlyReducedComparedToMaxTheta(max_theta_red_);
  filter_.utils_.out(NOX::Utils::Debug) << "Blocking::AddFilterIterate IsSufficientlyReduced = "
                                        << (is_sufficiently_red ? "TRUE" : "FALSE") << "\n";

  if (filter_iterates_.size() > 0)
  {
    auto& last = filter_iterates_.back();
    const int diff = newton_iter - last.first;
    if (diff == 0)
    {
      if (is_sufficiently_red) ++last.second;
      return;
    }

    if (diff != 1)
    {
      filter_iterates_.clear();
    }
  }

  if (is_sufficiently_red)
    filter_iterates_.push_back(std::make_pair(static_cast<unsigned>(newton_iter), 1));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Blocking::ReinitializeFilter()
{
  // get the number of rejections due to the filter in the last Newton iterate
  const unsigned last_rej_num = filter_iterates_.back().second;
  const bool do_reinit =
      (filter_iterates_.size() >= consecutive_iter_ or last_rej_num > consecutive_ls_steps_);

  if (do_reinit)
  {
    PrintInfo(filter_.utils_.out());

    Point::reinitFilter(filter_.filter_, filter_.theta_, max_theta_red_);
    filter_iterates_.clear();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Blocking::PrintInfo(std::ostream& os) const
{
  os << NOX::Utils::fill(40, '+') << "\n";
  os << "FILTER REINITIALIZATION"
     << "\n";
  os << NOX::Utils::fill(40, '+') << "\n";
  os << "The filter is going to be reinitialized due to the "
        "collected information over the last iterates.\n"
        "The filter rejected the following consecutive Newton iterates:\n\n";

  os << "{ Iteration, Number of Rejections }\n";
  for (const auto& p : filter_iterates_) os << "{ " << p.first << ", " << p.second << "}\n";

  os << NOX::Utils::fill(40, '+') << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::Filter::PostCheckStatus(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Solver::Generic& solver,
    const NOX::Abstract::Group& grp, NOX::StatusTest::CheckType checkType)
{
  const int iter_ls = linesearch.GetNumIterations();
  if (iter_ls == 1 and status_ == status_step_too_long)
  {
    return soc_->execute(linesearch, solver, grp, checkType);
  }
  else if (blocking_.filter_iterates_.size())
  {
    blocking_.ReinitializeFilter();
  }
  else
  {
    ThrowIfStepTooShort(linesearch, solver);
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::ThrowIfStepTooShort(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Solver::Generic& solver) const
{
  const double& step = linesearch.GetStepLength();
  const enum NOX::StatusTest::StatusType active_set_status = GetActiveSetStatus(solver);

  if (status_ == status_step_too_short)
    dserror(
        "The step-length %1.6e is lower than the minimal allowed step length"
        "amin = %1.6e. This indicates that we can't find a feasible solution "
        "in the current search direction and we decide to stop here. (active-set status = %s)",
        step, gamma_alpha_ * amin_,
        NOX::NLN::StatusTest::StatusType2String(active_set_status).c_str());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::SecondOrderCorrection::execute(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Solver::Generic& solver,
    const NOX::Abstract::Group& grp, NOX::StatusTest::CheckType checkType)
{
  // avoid recursive execution calls
  if (issoc_) return filter_.GetStatus();
  issoc_ = true;

  NOX::Abstract::Group& mutable_grp = const_cast<NOX::Abstract::Group&>(grp);
  NOX::NLN::Group& nln_grp = dynamic_cast<NOX::NLN::Group&>(mutable_grp);
  curr_type_ = whichType(solver);

  filter_.utils_.out() << "SOC-type       = " << CorrectionType2String(curr_type_) << "\n";

  {
    const double t_start = Teuchos::Time::wallTime();
    computeSystem(nln_grp, solver);
    solve(linesearch, solver, nln_grp);
    time_exe_ = Teuchos::Time::wallTime() - t_start;
  }

  postprocess(linesearch, solver, nln_grp, checkType);
  issoc_ = false;

  return soc_status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CorrectionType NOX::NLN::INNER::StatusTest::Filter::SecondOrderCorrection::whichType(
    const NOX::Solver::Generic& solver) const
{
  switch (user_type_)
  {
    // if an user-specified SOC type is provided we use that one
    case CorrectionType::soc_cheap:
    case CorrectionType::soc_full:
      return user_type_;
    // let the algorithm decide
    case CorrectionType::soc_automatic:
      return automaticTypeChoice(solver);
    // no SOC correction type
    default:
    {
      dserror(
          "Unsupported Second Order Correction type detected. "
          "Given type = \"%s\"",
          CorrectionType2String(user_type_).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CorrectionType
NOX::NLN::INNER::StatusTest::Filter::SecondOrderCorrection::automaticTypeChoice(
    const NOX::Solver::Generic& solver) const
{
  const enum NOX::StatusTest::StatusType active_set_status = filter_.GetActiveSetStatus(solver);
  const Point& trial_fp = *filter_.curr_fpoints_.second;
  const bool isfeasible = trial_fp.IsFeasible(filter_.GetConstraintTolerance(solver));

  switch (active_set_status)
  {
    /* If the active-set seems to be converged, we will perform the cheap
     * correction.
     *
     * EDIT: If the current filter point is feasible, i.e. the constraint
     * violation is below the given tolerance, a further reduction might be
     * hard to achieve and a full soc step is the better choice. The same is
     * true, if the bodies are out of contact, then the cheap step won't change
     * anything. */
    case NOX::StatusTest::Converged:
      return (isfeasible ? CorrectionType::soc_full : CorrectionType::soc_cheap);
    /* If the active-set is not yet converged or the status is undefined, a full
     * step will be performed. */
    default:
      return CorrectionType::soc_full;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::SecondOrderCorrection::postprocess(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Solver::Generic& solver,
    NOX::Abstract::Group& grp, NOX::StatusTest::CheckType checkType)
{
  soc_status_ = status_unevaluated;
  try
  {
    grp.computeF();
    // Note that this tests more than just the filter method. Also all
    // pre-testings will be executed.
    soc_status_ = linesearch.CheckInnerStatus(solver, grp, checkType);
  }
  catch (const char* e)
  {
    soc_status_ = status_step_too_long;
    filter_.utils_.out(NOX::Utils::Warning)
        << "WARNING: computeF failed "
           "after the SOC step. Recover and start line-search...\n";
    feclearexcept(FE_ALL_EXCEPT);
  }

  if (soc_status_ != status_converged)
  {
    const double t_start = Teuchos::Time::wallTime();
    filter_.RecoverFromBackup(grp);
    time_recover_ = Teuchos::Time::wallTime() - t_start;
  }

  print(filter_.utils_.out());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::SecondOrderCorrection::solve(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Solver::Generic& solver,
    NOX::Abstract::Group& grp) const
{
  // copy parameter list and perform the second order correction step
  Teuchos::ParameterList pnewton(
      solver.getList().sublist("Direction").sublist("Newton").sublist("Linear Solver"));

  // construct a new direction pointer
  const NOX::Abstract::Vector& x = grp.getX();
  const Epetra_BlockMap& map = dynamic_cast<const NOX::Epetra::Vector&>(x).getEpetraVector().Map();

  Teuchos::RCP<Epetra_Vector> dir_ptr = Teuchos::rcp(new Epetra_Vector(map, true));
  Teuchos::RCP<NOX::Epetra::Vector> nox_dir =
      Teuchos::rcp(new NOX::Epetra::Vector(dir_ptr, NOX::Epetra::Vector::CreateView));

  // compute the new direction
  const NOX::NLN::Solver::LineSearchBased& nln_solver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased&>(solver);
  NOX::Direction::Generic& direction = nln_solver.GetDirection();
  bool success = direction.compute(*nox_dir, grp, solver);
  if (not success) dserror("Solving of the SOC system failed!");

  // update the state in the group object
  grp.computeX(grp, *nox_dir, linesearch.GetStepLength());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::SecondOrderCorrection::computeSystem(
    NOX::NLN::Group& grp, const NOX::Solver::Generic& solver) const
{
  grp.computeCorrectionSystem(curr_type_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::SecondOrderCorrection::print(std::ostream& os) const
{
  os << "\n" << NOX::Utils::fill(46, '=') << "\n";
  os << "Performed a Second Order Correction (SOC) step\n";
  os << "User SOC-type  = " << CorrectionType2String(user_type_) << "\n";
  os << "SOC-type       = " << CorrectionType2String(curr_type_) << "\n";
  os << "Execution time = " << time_exe_ << " sec.\n";
  os << "New status     = " << StatusType2String(soc_status_) << "\n";
  if (soc_status_ != status_converged) os << "Recover time   = " << time_recover_ << " sec.\n";
  os << NOX::Utils::fill(46, '=') << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::BackupState::create(
    const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& dir)
{
  const NOX::NLN::Group* nln_grp_ptr = dynamic_cast<const NOX::NLN::Group*>(&grp);
  if (not nln_grp_ptr) dserror("Dynamic_cast failed!");

  if (xvector_.is_null())
  {
    const NOX::Epetra::Vector& nox_epetra_x = dynamic_cast<const NOX::Epetra::Vector&>(grp.getX());
    xvector_ = Teuchos::rcp(new NOX::Epetra::Vector(nox_epetra_x));
  }
  else
    xvector_->update(1.0, grp.getX(), 0.0);

  normf_ = grp.getF().norm(NOX::Abstract::Vector::TwoNorm);

  nln_grp_ptr->CreateBackupState(dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::BackupState::recover(NOX::Abstract::Group& grp) const
{
  if (xvector_.is_null()) dserror("The xvector_ is not set!");

  grp.setX(*xvector_);

  NOX::NLN::Group& nln_grp = dynamic_cast<NOX::NLN::Group&>(grp);
  nln_grp.RecoverFromBackupState();

  nln_grp.computeF();

  checkRecoveredState(nln_grp.getF());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::BackupState::checkRecoveredState(
    const NOX::Abstract::Vector& f) const
{
  const double recovered_normf = f.norm(NOX::Abstract::Vector::TwoNorm);
  static const double rtol = std::numeric_limits<double>::epsilon();

  const double ratiof = recovered_normf / normf_;

  if (ratiof < 1.0 - rtol or ratiof > 1.0 + rtol)
    dserror(
        "The recovery from the previously stored state failed! "
        "[recovery_normf / backup_normf = %.15e ]",
        ratiof);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::RecoverFromBackup(NOX::Abstract::Group& grp) const
{
  backup_.recover(grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::Filter::IsAdmissibleStep(
    const NOX::Solver::Generic& solver, const double& step) const
{
  if (step >= gamma_alpha_ * amin_)
    return status_step_too_long;
  else
    return status_step_too_short;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::StatusTest::StatusType NOX::NLN::INNER::StatusTest::Filter::GetActiveSetStatus(
    const NOX::Solver::Generic& solver) const
{
  const NOX::NLN::Solver::LineSearchBased* ls_solver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&solver);
  if (not ls_solver) dserror("The given non-linear solver is not line search based!");

  NOX::StatusTest::Generic* active_set_test =
      ls_solver->GetOuterStatusTest<NOX::NLN::StatusTest::ActiveSet>();

  if (not active_set_test) return NOX::StatusTest::Unevaluated;

  return active_set_test->checkStatus(solver, NOX::StatusTest::Complete);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::SufficientReductionCheck(const Point& trial_fp) const
{
  const Point& previous_fp = *curr_fpoints_.first;

  for (unsigned i = 0; i < Point::num_coords_; ++i)
  {
    if (trial_fp(i) <= previous_fp(i) - previous_fp.margin_(i))
    {
      return status_converged;
    }
  }

  return status_step_too_long;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::AugmentFilter()
{
  const Teuchos::RCP<Point> new_fp_ptr = curr_fpoints_.second;

  std::fill(filter_.begin(), filter_.end(), Teuchos::null);
  filter_.resize(non_dominated_filter_points_.size() + 1, Teuchos::null);
  plain_point_set::iterator it = filter_.begin();

  if (non_dominated_filter_points_.empty())
  {
    filter_[0] = new_fp_ptr;
  }
  else
  {
    bool is_inserted = false;
    for (plain_point_set::const_iterator cit = non_dominated_filter_points_.begin();
         cit != non_dominated_filter_points_.end(); ++cit, ++it)
    {
      const Teuchos::RCP<Point>& fp_ptr = *cit;

      if (!is_inserted and new_fp_ptr->norm_ < fp_ptr->norm_)
      {
        *(it++) = new_fp_ptr;
        is_inserted = true;
      }

      *it = *cit;
    }

    if (!is_inserted) filter_.back() = new_fp_ptr;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::FilterStatusType
NOX::NLN::INNER::StatusTest::Filter::AcceptabilityCheck(const Point& trial_fp)
{
  const unsigned prefiltering_index = Prefiltering(trial_fp);

  /* If prefiltering did not succeed, we perform the acceptability check point
   * by point */
  // if the current filter is empty the check passes by default
  bool passed_check = true;
  for (auto cit = filter_.begin(); cit != filter_.begin() + prefiltering_index; ++cit)
  {
    const Point& fp = **cit;

    passed_check = false;
    for (unsigned i = 0; i < fp.num_coords_; ++i)
    {
      if (fp.norm_ < 0.0 or trial_fp(i) < fp(i) - fp.margin_(i))
      {
        passed_check = true;
        break;
      }
    }

    // If the check failed for one filter point, the whole test failed.
    if (not passed_check)
    {
      utils_.out(NOX::Utils::Debug) << "\n" << NOX::Utils::fill(22, '=') << "\n";
      utils_.out(NOX::Utils::Debug) << "rejected by:\n";
      fp.print(utils_.out(NOX::Utils::Debug));
      utils_.out(NOX::Utils::Debug) << "\n" << NOX::Utils::fill(22, '=') << "\n";
      break;
    }
  }

  return (passed_check ? FilterStatusType::passed_point_by_point : FilterStatusType::rejected);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned NOX::NLN::INNER::StatusTest::Filter::Prefiltering(const Point& trial_fp)
{
  const double sqrt_num_coords =
      std::sqrt<double>(static_cast<double>(trial_fp.num_coords_)).real();
  const double max_gamma = std::max(Point::gamma_obj_, Point::gamma_theta_);

  /* The following pre-filtering approach does only work, if ALL involved
   * coordinates are positive (or equal to zero). Therefore, if one of the trial
   * filter coordinates is negative, we will skip this pre-filtering from
   * thereon. */
  static bool skip_prefiltering = false;
  const double* min_coord =
      std::min_element(trial_fp.coords_.A(), trial_fp.coords_.A() + trial_fp.num_coords_);
  if (!skip_prefiltering and *min_coord < 0.0)
  {
    skip_prefiltering = true;
    utils_.out(NOX::Utils::Debug) << "\nA negative coordinate has been detected "
                                     "[val = "
                                  << NOX::Utils::sciformat(*min_coord, OUTPUT_PRECISION)
                                  << "]. The pre-filtering will be skipped from "
                                     "now on...\n\n";
  }

  unsigned prefiltering_index = 0;
  unsigned non_dominated_index = 0;

  if (skip_prefiltering)
  {
    prefiltering_index = filter_.size();
    non_dominated_index = 0;
  }
  else
  {
    for (plain_point_set::const_iterator cit = filter_.begin(); cit < filter_.end();
         ++cit, ++prefiltering_index)
    {
      const Point& fp = **cit;

      if (fp.norm_ > 0.0 and trial_fp.norm_ > fp.norm_) ++non_dominated_index;

      if (fp.norm_ < 0.0 or trial_fp.norm_ < fp.norm_ - sqrt_num_coords * max_gamma * fp.maxTheta())
        break;
    }
  }

  non_dominated_filter_points_.clear();
  IdentifyNonDominatedFilterPoints(trial_fp, non_dominated_index);

  return prefiltering_index;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::IdentifyNonDominatedFilterPoints(
    const Point& trial_fp, const unsigned non_dominated_index)
{
  non_dominated_filter_points_.reserve(filter_.size());

  std::copy(filter_.begin(), filter_.begin() + non_dominated_index,
      std::back_inserter(non_dominated_filter_points_));

  for (plain_point_set::const_iterator cit = filter_.begin() + non_dominated_index;
       cit < filter_.end(); ++cit)
  {
    const Point& fp = **cit;

    unsigned num_dominated_coords = 0;

    if (fp.norm_ < 0.0)
      num_dominated_coords = fp.num_coords_;
    else
    {
      for (unsigned i = 0; i < fp.num_coords_; ++i)
      {
        if (fp(i) >= trial_fp(i)) ++num_dominated_coords;
      }
    }

    // If the trial filter point does not dominate the current filter point, we
    // will keep the current filter point in our filter point set.
    if (num_dominated_coords < fp.num_coords_)
    {
      non_dominated_filter_points_.push_back(*cit);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::SetupModelTerms(const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp, const Interface::Required& interface)
{
  const NOX::MeritFunction::Generic& merit_func = interface.GetMeritFunction();
  if (dynamic_cast<const MeritFunction::Lagrangian*>(&merit_func))
  {
    const NOX::NLN::MeritFunction::Lagrangian& lagrangian =
        dynamic_cast<const NOX::NLN::MeritFunction::Lagrangian&>(merit_func);

    model_lin_terms_(0) = lagrangian.computeSlope(dir, grp);
    model_mixed_terms_(0) = lagrangian.computeMixed2ndOrderTerms(dir, grp);
  }
  else
    dserror("Currently unsupported merit function type: \"%s\"", merit_func.name().c_str());

  theta_.computeSlope(dir, grp, model_lin_terms_.A() + Point::num_obj_coords_);
  theta_.computeMixed2ndOrderTerms(dir, grp, model_mixed_terms_.A() + Point::num_obj_coords_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::ComputeMinimalStepLengthEstimates()
{
  /* compute minimal step length estimate based on the 2nd objective function
   * filter acceptability check */
  amin_obj_ = MinimalStepLengthEstimateOfObjFuncFilterCheck();

  /* compute minimal step length estimate based on the 2nd constraint violation
   * filter acceptability check */
  amin_theta_ = theta_.minimalStepLengthEstimate(curr_points_.first->A() + Point::num_obj_coords_,
      model_lin_terms_.A() + Point::num_obj_coords_);

  /* compute minimal step length estimate based on the f-type switching condition */
  amin_ftype_ = 1.0;
  if (CheckFTypeSwitchingCondition(1.0))
  {
    amin_ftype_ = MinimalStepLengthEstimateOfFTypeCondition();
  }

  amin_ = std::min(std::min(amin_obj_, amin_theta_), amin_ftype_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::MinimalStepLengthEstimateOfObjFuncFilterCheck() const
{
  double amin_obj = 1.0;

  // Is the current search direction a descent direction for the objective model?
  if (model_lin_terms_(0) < 0.0)
  {
    // get the maximal value of the accepted infeasibility measurements
    const double max_theta = curr_points_.first->maxTheta();
    // get the accepted objective function value
    const double obj_slope = model_lin_terms_(0);

    // check the 2nd order mixed derivative term
    const double obj_mixed_term = model_mixed_terms_(0);
    const bool is_linear_obj_model = (std::abs(obj_mixed_term) < 1.0e-12);

    if (is_linear_obj_model)
    {
      /*-----------------------------*
       | Filter Check (linear model) |
       *-----------------------------*-------------------------------------*
       |  Linear model:                                                    |
       |  s_f* (L_k + a * LIN(L_k)) < s_f * L_k - s_t * gamma_f * theta_k, |
       |                                                                   |
       |  where s_f and s_t are the scaling factors for the objective      |
       |  and constraint values, respectively. Note that we assume that    |
          LIN(L_k) is negative, i.e. descent direction.                    |
       |                                                                   |
       | => a > - (s_t/s_f) * (gamma_f * theta_k)/LIN(L_k).                |
       *-------------------------------------------------------------------*/

      amin_obj =
          -(Point::gamma_obj_ * max_theta * Point::scale_(1)) / (obj_slope * Point::scale_(0));
    }
    else
    {
      /*-------------------------------*
       | Filter Check (2nd order model |
       *-------------------------------*---------------------------------------------*
       |  Quadratic model:                                                           |
       |  s_f * (L_k + c1 * a + c2 * a^2) < s_f * L_k - s_t * gamma_f * theta_k      |
       |                                                                             |
       |  a_1/2 =   (-c1 (+-) sqrt(c1^2 + 4 * c2 * gamma_f * (s_t/s_f) * theta_k))   |
       |          / (-2 * c2).                                                       |
       |  Only the solution corresponding to the minus sign is interesting. To       |
       |  understand this, we consider two different cases. For all of them is       |
       |  c1 lower than zero (descent direction):                                    |
       |                                                                             |
       |  [1] c2 > 0. This corresponds to a parabola which opens upward:             |
       |      In this case there are normally two positive roots and we choose the   |
       |      1st/smaller one. The minimizer of the quadratic 1-D model is not       |
       |      important for us. Nevertheless, it would be possible to check the      |
       |      gradient of the 1-D model for the unity step length and extend the     |
       |      line search method by increasing the step-length if the gradient       |
       |                                                                             |
       |  r_s(x_k+d) - (z_n + dz)^T * grad[wgn(x_k+d)]^T * d - wgn(x_k+d)^T * dz     |
       |                                                                             |
       |      is lower than zero and the step is not accepted. At the moment we use  |
       |      always a backtracking strategy and need only a lower bound for the     |
       |      step length parameter.                                                 |
       |                                                                             |
       |  [2] c2 < 0. This corresponds to a parabola which opens downward. Here we   |
       |      use the same idea. We are interested in a lower bound. We need the     |
       |      2nd/right root, which corresponds to the minus sign again.             |
       *-----------------------------------------------------------------------------*/
      // If the parabola opens upward and the minimum of the quadratic model lies over
      // the specified threshold, it is not possible to find a solution, because the
      // the parabola and the constant line have no intersection point.
      if (obj_mixed_term > 0.0 and
          ((obj_slope * obj_slope) / (-4.0 * obj_mixed_term) >
              -Point::scale_(1) / Point::scale_(0) * Point::gamma_obj_ * max_theta))
        amin_obj = 1.0;
      else
      {
        double atmin_obj = std::sqrt<double>(
            obj_slope * obj_slope - 4.0 * obj_mixed_term * Point::gamma_obj_ * Point::scale_(1) /
                                        Point::scale_(0) * max_theta)
                               .real();
        atmin_obj = (-obj_slope - atmin_obj) / (2.0 * obj_mixed_term);

        amin_obj = std::min(amin_obj, std::max(0.0, atmin_obj));
      }
    }
  }

  return amin_obj;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::MinimalStepLengthEstimateOfFTypeCondition() const
{
  // linear term of the objective model
  const double obj_slope = model_lin_terms_(0);

  // 2nd order mixed derivative term of the objective function
  const double obj_mixed = model_mixed_terms_(0);

  const double scale_of_max_theta = curr_points_.first->scaleOfMaxTheta();

  // safe-guarding strategy: lower and upper bounds for the minimal step length estimate
  double lBound = 0.0;
  double uBound = 1.0;

  const double d = std::pow(scale_of_max_theta, st_) / std::pow(Point::scale_(0), sf_);

  const double f_lbound = computeFTypeSwitchingCondition(lBound, d);
  if (f_lbound > 0.0) dserror("The function value for the lower bound is greater than zero!");

  const double f_ubound = computeFTypeSwitchingCondition(uBound, d);
  if (f_ubound < 0.0) dserror("The function value for the upper bound is lower than zero!");

  // set initial value (secant value)
  double amin = f_lbound / (f_lbound - f_ubound);

  // newton control parameters
  static const double TOL_LOCAL_NEWTON = 1.0e-8;
  static const unsigned ITERMAX = 50;
  unsigned iter = 0;
  bool isconverged = false;

  while (not isconverged and iter < ITERMAX)
  {
    const double f = computeFTypeSwitchingCondition(amin, d);

    // update lower bound
    if (f < 0.0 and amin > lBound) lBound = amin;
    // update upper bound
    else if (f > 0.0 and amin < uBound)
      uBound = amin;

    double da = -(obj_slope + obj_mixed * amin);
    da = -f / (std::pow(da, sf_) * (1.0 - sf_ * obj_mixed * amin / (da)));

    amin += da;

    // safe-guarding strategy
    if (amin < lBound or amin > uBound) amin = 0.5 * (lBound + uBound);

    // relative convergence check
    isconverged = std::abs(da) < TOL_LOCAL_NEWTON * std::max(amin, 1.0e-12);

    ++iter;
  }
  if (not isconverged) dserror("The local Newton did not converge in %d steps!", iter);

  return amin;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::INNER::StatusTest::Filter::CheckFTypeSwitchingCondition(const double step) const
{
  const double obj_model = GetObjModel(step);

  // Is the previously accepted infeasibility measure small enough and the
  // current direction a descent direction for the used objective function?
  if ((curr_points_.first->maxTheta() < theta_min_ftype_) and (obj_model < 0.0))
  {
    const double scale_of_max_theta = curr_points_.first->scaleOfMaxTheta();

    const double d = std::pow(scale_of_max_theta, st_) / std::pow(Point::scale_(0), sf_);

    return (computeFTypeSwitchingCondition(step, d) > 0.0);
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::computeFTypeSwitchingCondition(
    const double step, const double d) const
{
  dsassert(d > 0.0, "The scaling factor d is smaller than / equal to zero!");

  const double max_theta = curr_points_.first->maxTheta();

  // linear term of the objective model
  const double obj_slope = model_lin_terms_(0);

  // 2nd order mixed derivative term of the objective function
  const double obj_mixed = model_mixed_terms_(0);

  utils_.out() << "obj_slope   = " << obj_slope << "\n"
               << "obj_mixed   = " << obj_mixed << "\n"
               << "step        = " << step << "\n"
               << "d           = " << d << "\n"
               << "max_theta   = " << max_theta << "\n"
               << "theta_min_f = " << theta_min_ftype_ << "\n";

  utils_.out() << "f-type-result = "
               << std::pow(-(obj_slope + step * obj_mixed), sf_) * step -
                      d * std::pow(max_theta, st_)
               << "\n";

  return std::pow(-(obj_slope + step * obj_mixed), sf_) * step - d * std::pow(max_theta, st_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::GetObjModel(const double step) const
{
  return step * model_lin_terms_(0) + step * step * model_mixed_terms_(0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::Infeasibility::minimalStepLengthEstimate(
    const double* accepted_theta, const double* theta_slope) const
{
  double amin = 1.0;
  /*------------------*
   | Filter Check     |
   *------------------*------------------------------------*
   |    theta_k + a * LIN(theta_k) < (1-gamma_t) * theta_k |
   |                                                       |
   | => a > -gamma_t*theta_k / LIN(theta_k)                |
   *-------------------------------------------------------*/
  for (unsigned i = 0; i < number_; ++i)
  {
    double amin_theta = 1.0;

    /* Is the current search direction a descent direction for the infeasibility
     * measure? */
    if (theta_slope[i] < 0.0)
    {
      amin_theta = -Point::gamma_theta_ * accepted_theta[i] / theta_slope[i];
    }

    amin = std::min(amin_theta, amin);
  }

  return amin;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Infeasibility::computeSlope(
    const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp,
    double* theta_slope_values) const
{
  plain_merit_func_set::const_iterator cit = vector_.begin();

  for (unsigned i = 0; cit != vector_.end(); ++cit, ++i)
  {
    const NOX::MeritFunction::Generic& func = **cit;

    theta_slope_values[i] = func.computeSlope(dir, grp);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Infeasibility::computeMixed2ndOrderTerms(
    const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp,
    double* theta_mixed_values) const
{
  // no mixed 2nd order terms for the infeasibility measures
  std::fill(theta_mixed_values, theta_mixed_values + number_, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::Filter::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::INNER::StatusTest::Filter::Print(std::ostream& stream, int indent) const
{
  std::string indent_str;
  indent_str.assign(indent, ' ');

  stream << indent_str;
  stream << status_;

  std::string par_indent("    ");
  par_indent += indent_str;
  const int par_length = par_indent.size();

  stream << indent_str << "The filter " << FilterStatus2String(filter_status_) << ".\n";
  // Skip the remaining output if the filter check has not been evaluated.
  if (status_ == status_unevaluated) return stream;

  stream << par_indent << "CURRENT POINT PAIR \n{\n";
  stream << par_indent << "+ Previously accepted point:\n";
  curr_points_.first->print(stream, par_length, utils_);

  stream << par_indent << "+ Current trial point:\n";
  curr_points_.second->print(stream, par_length, utils_);
  stream << "}\n";

  stream << indent_str << "CURRENT FILTER POINT PAIR \n{\n";
  stream << par_indent << "+ Previously accepted filter point:\n";
  curr_fpoints_.first->print(stream, par_length, utils_);

  stream << par_indent << "+ Current trial filter point:\n";
  curr_fpoints_.second->print(stream, par_length, utils_);
  stream << "}\n";

  stream << indent_str << "FILTER\n{\n";
  unsigned id = 0;
  for (auto cit = filter_.begin(); cit != filter_.end(); ++cit, ++id)
  {
    const Point& fp = **cit;
    stream << "(" << id << ") ";
    fp.print(stream, par_length, utils_);
  }
  stream << "}\n";

  stream << indent_str << "F-Type condition is " << (is_ftype_step_ ? "" : "not ")
         << "fulfilled.\n";
  if (is_ftype_step_) armijo_test_->Print(stream, par_indent.size());

  if (not utils_.isPrintType(NOX::Utils::Details)) return stream;

  // --- Detailed filter output -----------------------------------------------

  stream << indent_str << "MINIMAL STEP LENGTH ESTIMATES\n{\n";
  stream << par_indent
         << "Objective estimate = " << NOX::Utils::sciformat(amin_obj_, OUTPUT_PRECISION) << "\n";
  stream << par_indent
         << "Theta estimate     = " << NOX::Utils::sciformat(amin_theta_, OUTPUT_PRECISION) << "\n";
  stream << par_indent
         << "F-type estimate    = " << NOX::Utils::sciformat(amin_ftype_, OUTPUT_PRECISION) << "\n";
  stream << par_indent << "-------------------- "
         << "\n";
  stream << par_indent << "Over-all estimate  = " << NOX::Utils::sciformat(amin_, OUTPUT_PRECISION)
         << "\n";
  stream << "}\n";

  stream << indent_str << "INFEASIBILITY STATISTICS\n{\n";
  stream << par_indent << "Number of theta  = " << theta_.number_ << "\n";
  stream << par_indent << "Types            = {";
  for (const auto& theta_ptr : theta_.vector_)
  {
    stream << " \"" << theta_ptr->name() << "\"";
  }
  stream << " };\n";
  stream << "}\n";

  stream << indent_str << "GENERAL POINT STATISTICS\n{\n";
  stream << par_indent << "Number of coords = " << Point::num_coords_ << "\n";
  stream << par_indent << "Number of obj    = " << Point::num_obj_coords_ << "\n";
  stream << par_indent
         << "Gamma_obj        = " << NOX::Utils::sciformat(Point::gamma_obj_, OUTPUT_PRECISION)
         << "\n";
  stream << par_indent
         << "Gamma_theta      = " << NOX::Utils::sciformat(Point::gamma_theta_, OUTPUT_PRECISION)
         << "\n";
  stream << par_indent << "Scales           = { ";
  for (unsigned i = 0; i < Point::num_coords_; ++i)
  {
    if (i != 0) stream << ", ";
    stream << NOX::Utils::sciformat(Point::scale_(i), OUTPUT_PRECISION);
  }
  stream << " };\n";
  stream << par_indent << "Valid scales     = {";
  for (const bool valid_scale : Point::isvalid_scaling_)
  {
    stream << " ";
    stream << (valid_scale ? "valid" : "invalid");
  }
  stream << " };\n";
  stream << par_indent << "Weights          = { ";
  for (unsigned i = 0; i < Point::num_coords_; ++i)
  {
    if (i != 0) stream << ", ";
    stream << NOX::Utils::sciformat(Point::weights_(i), OUTPUT_PRECISION);
  }
  stream << " };\n";
  stream << par_indent << "Global Max Theta = { ";
  for (int i = 0; i < Point::global_scaled_max_thetas_.Length(); ++i)
  {
    if (i != 0) stream << ", ";
    stream << NOX::Utils::sciformat(Point::global_scaled_max_thetas_(i), OUTPUT_PRECISION);
  }
  stream << " }; [already scaled]\n";
  stream << par_indent << "Initial Global Max Theta Scale = " << Point::global_init_max_theta_scale_
         << "\n";
  stream << "}\n";

  stream << indent_str << "FILTER BLOCKING STATISTICS\n{\n";
  stream << par_indent << "Number of entries   = " << blocking_.filter_iterates_.size() << "\n";
  stream << par_indent << "Entries ( IT; NUM ) = {";
  for (const auto& p : blocking_.filter_iterates_)
  {
    stream << " { " << p.first << ", " << p.second << " } ";
  }
  stream << " };\n";
  stream << "}\n";

  return stream;
}
