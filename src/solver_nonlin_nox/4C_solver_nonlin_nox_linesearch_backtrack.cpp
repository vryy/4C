// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_linesearch_backtrack.hpp"  // class definition

#include "4C_comm_mpi_utils.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_linesearch_prepostoperator.hpp"
#include "4C_solver_nonlin_nox_solver_linesearchbased.hpp"
#include "4C_solver_nonlin_nox_statustest_normf.hpp"
#include "4C_solver_nonlin_nox_vector.hpp"
#include "4C_utils_exceptions.hpp"

#include <fenv.h>
#include <mpi.h>
#include <NOX_GlobalData.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef FOUR_C_ENABLE_FE_TRAPPING
#include <cfenv>
#endif

FOUR_C_NAMESPACE_OPEN

namespace
{



  /**
   * @brief Runs the \p evaluation and detects whether floating point exceptions occurred on any
   * rank.
   *
   * @param evaluation The evaluation to run.
   * @param allow_floating_point_exceptions Whether to allow floating point exceptions or not.
   * @param comm The MPI communicator to use for the reduction. No reduction is performed if \c
   * FOUR_C_ENABLE_FE_TRAPPING is not set or if \p allow_floating_point_exceptions is false.
   * @param os The output stream to use for logging rank-local floating point exceptions.
   * @return true if a floating point exception occurred on any rank, false otherwise.
   */
  bool run_and_detect_floating_point_exceptions(const std::function<void()>& evaluation,
      bool allow_floating_point_exceptions, MPI_Comm comm, std::ostream& os)
  {
#ifdef FOUR_C_ENABLE_FE_TRAPPING
    if (allow_floating_point_exceptions)
    {
      // Floating point exceptions for which backtracking should be performed.
      constexpr int handled_floating_point_exceptions = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW;

      /// Guard to set the floating point exception environment and restore it after evaluation.
      class FloatingPointExceptionGuard
      {
       public:
        explicit FloatingPointExceptionGuard()
        {
          fedisableexcept(handled_floating_point_exceptions);
          // clear any left-over flags
          feclearexcept(handled_floating_point_exceptions);
        }

        FloatingPointExceptionGuard(const FloatingPointExceptionGuard&) = delete;
        FloatingPointExceptionGuard& operator=(const FloatingPointExceptionGuard&) = delete;

        ~FloatingPointExceptionGuard()
        {
          feclearexcept(handled_floating_point_exceptions);
          // restore the previously enabled floating point exceptions
          fedisableexcept(FE_ALL_EXCEPT);
          feenableexcept(previously_enabled_floating_point_exceptions_);
        }

       private:
        const int previously_enabled_floating_point_exceptions_ = fegetexcept();
      };

      FloatingPointExceptionGuard guard;

      evaluation();

      const bool local_fe_except = fetestexcept(handled_floating_point_exceptions) != 0;
      feclearexcept(handled_floating_point_exceptions);

      if (local_fe_except)
      {
        os << "WARNING: Floating point exception occurred on rank "
           << Core::Communication::my_mpi_rank(comm);
      }

      return Core::Communication::sum_all(static_cast<int>(local_fe_except), comm) > 0;
    }
#endif

    evaluation();
    return false;
  }
}  // namespace

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::LineSearch::Backtrack::Backtrack(const Teuchos::RCP<::NOX::GlobalData>& gd,
    const Teuchos::RCP<::NOX::StatusTest::Generic> outerTests,
    const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic> innerTests,
    Teuchos::ParameterList& params)
    : ls_iters_(0),
      step_ptr_(nullptr),
      default_step_(0.0),
      reduction_factor_(0.0),
      check_type_(::NOX::StatusTest::Complete),
      status_(NOX::Nln::Inner::StatusTest::status_unevaluated),
      outer_tests_ptr_(outerTests),
      inner_tests_ptr_(innerTests),
      pre_post_operator_ptr_(Teuchos::null),
      controller_(*this)
{
  reset(gd, params);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::LineSearch::Backtrack::reset(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  Teuchos::ParameterList& p = params.sublist("Backtrack");

  utils_ = gd->getUtils();
  merit_function_ptr_ = gd->getMeritFunction();

  ls_iters_ = 0;
  search_direction_ptr_ = Teuchos::null;

  status_ = NOX::Nln::Inner::StatusTest::status_unevaluated;

  default_step_ = p.get("Default Step", 1.0);
  reduction_factor_ = p.get("Reduction Factor", 0.5);
  if ((reduction_factor_ <= 0.0) || (reduction_factor_ >= 1.0))
  {
    std::ostringstream msg;
    msg << "Invalid choice \"" << reduction_factor_ << "\" for \"Reduction Factor\"!\n"
        << "Value must be greater than zero and less than 1.0.";
    throw_error("reset", msg.str());
  }

  check_type_ = Teuchos::getIntegralValue<::NOX::StatusTest::CheckType>(
      params, "Inner Status Test Check Type");

  allow_floating_point_exceptions_ = p.get("Allow Exceptions", false);

  pre_post_operator_ptr_ = Teuchos::make_rcp<PrePostOperator>(params);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::LineSearch::Backtrack::reset()
{
  ls_iters_ = 0;
  search_direction_ptr_ = Teuchos::null;

  status_ = NOX::Nln::Inner::StatusTest::status_unevaluated;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::LineSearch::Backtrack::compute(::NOX::Abstract::Group& grp, double& step,
    const ::NOX::Abstract::Vector& dir, const ::NOX::Solver::Generic& s)
{
  // find out the current communicator
  const auto comm = dynamic_cast<const NOX::Nln::Vector&>(dir).get_linalg_vector().get_comm();
  // -------------------------------------------------
  // (re)set important line search parameters
  // -------------------------------------------------
  reset();
  // get the old solution group
  const ::NOX::Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  // update the search direction pointer
  search_direction_ptr_ = Teuchos::rcpFromRef(dir);
  // set the step pointer to the inserted step variable
  step_ptr_ = &step;
  // reset the step length
  step = default_step_;
  // initialize the inner status test
  // ----------------------------------------------------------------------
  // BE CAREFUL HERE:
  // During the copy operation in Solver::LineSearchBased::step()
  // the current grp loses the ownership of the sharedLinearSystem. If we
  // want to access the jacobian, we have to use the oldGrp
  // (target of the copy process), instead.                hiermeier 08/15
  // ----------------------------------------------------------------------
  if (not oldGrp.isJacobian()) throw_error("compute()", "Ownership changed unexpectedly!");

  /* Setup the inner status test */
  status_ = inner_tests_ptr_->check_status(controller_, oldGrp, check_type_);

  // increase iteration counter after initialization
  ++ls_iters_;

  // -------------------------------------------------
  // update the solution vector and get a trial point
  // -------------------------------------------------
  grp.computeX(oldGrp, dir, step);
  ::NOX::Abstract::Group::ReturnType rtype = ::NOX::Abstract::Group::Ok;
  bool fpe_occurred = false;

  fpe_occurred = run_and_detect_floating_point_exceptions([&]() { rtype = grp.computeF(); },
      allow_floating_point_exceptions_, comm, utils_->out(::NOX::Utils::Warning));
  if (rtype != ::NOX::Abstract::Group::Ok) throw_error("compute", "Unable to compute F!");

  if (fpe_occurred)
  {
    utils_->out() << "Last step caused a floating point exception. Reducing step length...\n";
    status_ = NOX::Nln::Inner::StatusTest::status_step_too_long;
  }
  else
  {
    /* Safe-guarding of the inner status test:
     * If the outer NormF test is converged for a full step length,
     * we don't have to reduce the step length any further.
     * This additional check becomes necessary, because of cancellation
     * errors and related numerical artifacts. */
    // check the outer status test for the full step length
    outer_tests_ptr_->checkStatus(s, check_type_);

    const NOX::Nln::Solver::LineSearchBased& lsSolver =
        static_cast<const NOX::Nln::Solver::LineSearchBased&>(s);

    const ::NOX::StatusTest::StatusType ostatus =
        lsSolver.get_status<NOX::Nln::StatusTest::NormF>();

    /* Skip the inner status test, if the outer NormF test is
     * already converged! */
    if (ostatus == ::NOX::StatusTest::Converged)
    {
      return true;
    }
  }

  // -------------------------------------------------
  // print header if desired
  // -------------------------------------------------

  utils_->out(::NOX::Utils::InnerIteration) << "\n"
                                            << ::NOX::Utils::fill(72, '=') << "\n"
                                            << "-- Backtrack Line Search -- \n";

  if (not fpe_occurred)
  {
    status_ = inner_tests_ptr_->check_status(controller_, grp, check_type_);
    print_update(utils_->out(::NOX::Utils::InnerIteration));
  }
  // -------------------------------------------------
  // inner backtracking loop
  // -------------------------------------------------
  while (status_ == NOX::Nln::Inner::StatusTest::status_step_too_long)
  {
    // -------------------------------------------------
    // reduce step length
    // -------------------------------------------------
    pre_post_operator_ptr_->run_pre_modify_step_length(s, *this);
    step *= reduction_factor_;

    // -------------------------------------------------
    // - update the solution vector and get a trial point
    // - increase line search step counter
    // -------------------------------------------------
    grp.computeX(oldGrp, dir, step);
    ++ls_iters_;

    fpe_occurred = run_and_detect_floating_point_exceptions([&]() { rtype = grp.computeF(); },
        allow_floating_point_exceptions_, comm, utils_->out(::NOX::Utils::Warning));
    if (rtype != ::NOX::Abstract::Group::Ok) throw_error("compute", "Unable to compute F!");
    if (fpe_occurred)
    {
      utils_->out() << "Last step caused a floating point exception. Reducing step length...\n";
      status_ = NOX::Nln::Inner::StatusTest::status_step_too_long;
    }
    else
    {
      status_ = inner_tests_ptr_->check_status(controller_, grp, check_type_);
      print_update(utils_->out(::NOX::Utils::InnerIteration));
    }
  }
  // -------------------------------------------------
  // print footer if desired
  // -------------------------------------------------
  utils_->out(::NOX::Utils::InnerIteration) << ::NOX::Utils::fill(72, '=') << "\n";

  if (status_ == NOX::Nln::Inner::StatusTest::status_step_too_short)
    throw_error("compute()",
        "The current step is too short and no "
        "restoration phase is implemented!");
  else if (status_ == NOX::Nln::Inner::StatusTest::status_no_descent_direction)
    throw_error("compute()", "The given search direction is no descent direction!");

  return (status_ == NOX::Nln::Inner::StatusTest::status_converged ? true : false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::Nln::LineSearch::Backtrack::print_update(std::ostream& os) const
{
  // Print the status test parameters at each iteration if requested
  if (status_ == NOX::Nln::Inner::StatusTest::status_step_too_long)
  {
    os << ::NOX::Utils::fill(72, '-') << "\n";
    os << "-- Inner Status Test Results --\n";
    inner_tests_ptr_->print(os);
    os << ::NOX::Utils::fill(72, '-') << "\n";
  }
  // Print the final parameter values of the status test
  if (status_ != NOX::Nln::Inner::StatusTest::status_step_too_long)
  {
    os << ::NOX::Utils::fill(72, '-') << "\n";
    os << "-- Final Inner Status Test Results --\n";
    inner_tests_ptr_->print(os);
    os << ::NOX::Utils::fill(72, '-') << "\n";
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::Nln::LineSearch::Backtrack::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::Nln::LineSearch::Backtrack::" << functionName << " - " << errorMsg
      << std::endl;
  FOUR_C_THROW("{}", msg.str());
}

NOX::Nln::LineSearch::Backtrack::Controller::Controller(
    NOX::Nln::LineSearch::Backtrack& backtrack_solver)
    : backtrack_solver_(backtrack_solver)
{
}

int NOX::Nln::LineSearch::Backtrack::Controller::get_num_iterations() const
{
  return backtrack_solver_.ls_iters_;
}

const ::NOX::MeritFunction::Generic&
NOX::Nln::LineSearch::Backtrack::Controller::get_merit_function() const
{
  if (backtrack_solver_.merit_function_ptr_.is_null())
    backtrack_solver_.throw_error(
        "Controller::get_merit_function", "The merit function pointer is not initialized!");

  return *backtrack_solver_.merit_function_ptr_;
}

const ::NOX::Abstract::Vector& NOX::Nln::LineSearch::Backtrack::Controller::get_search_direction()
    const
{
  if (backtrack_solver_.search_direction_ptr_.is_null())
    backtrack_solver_.throw_error(
        "Controller::get_search_direction", "The search direction ptr is not initialized!");

  return *backtrack_solver_.search_direction_ptr_;
}

double NOX::Nln::LineSearch::Backtrack::Controller::get_step_length() const
{
  if (backtrack_solver_.step_ptr_ == nullptr)
    backtrack_solver_.throw_error("Controller::get_step_length", "Step pointer is nullptr!");

  return *backtrack_solver_.step_ptr_;
}

void NOX::Nln::LineSearch::Backtrack::Controller::set_step_length(double step)
{
  if (backtrack_solver_.step_ptr_ == nullptr)
    backtrack_solver_.throw_error("Controller::set_step_length", "Step pointer is nullptr!");

  *backtrack_solver_.step_ptr_ = step;
}

FOUR_C_NAMESPACE_CLOSE
