// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_inner_statustest_armijo.hpp"  // class definition

#include "4C_linalg_vector.hpp"
#include "4C_solver_nonlin_nox_linesearch_controller.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Abstract_Group.H>
#include <NOX_MeritFunction_Generic.H>
#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::Armijo::Armijo(
    const double& c_1, const bool& isMonotone, const std::size_t& maxHistSize)
    : status_(status_unevaluated),
      c_1_(c_1),
      fref_(0.0),
      fcurr_(0.0),
      slope_(0.0),
      step_(1.0),
      isMonotone_(isMonotone),
      maxHistSize_(maxHistSize),
      histVector_(std::deque<double>(0))
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Inner::StatusTest::Armijo::setup(
    const NOX::Nln::LineSearch::Controller& ls_controller, const ::NOX::Abstract::Group& grp)
{
  const ::NOX::MeritFunction::Generic& mrtFct = ls_controller.get_merit_function();

  // get the reference merit function value
  fref_ = mrtFct.computef(grp);

  // get the slope once (doesn't change during the inner iteration)
  slope_ = mrtFct.computeSlope(ls_controller.get_search_direction(), grp);

  // return false if the search direction is no descent direction
  if (slope_ >= 0.0) return false;

  // -------------------------------------------
  // Non-monotone setup
  // -------------------------------------------
  if (not isMonotone_)
  {
    // Compare the current merit function value with the last
    // added value of the history vector and augment the history
    // vector only if the two values differ in more than
    // machine precision.
    if (fref_ != histVector_.front())
    {
      histVector_.push_front(fref_);
      // remove the last element if the maximum size of the history vector
      // is reached
      if (histVector_.size() > maxHistSize_) histVector_.pop_back();
    }
    // get the maximal merit function value of the last accepted steps for
    // the Armijo check
    fref_ = *std::max_element(histVector_.begin(), histVector_.end());
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::StatusType NOX::Nln::Inner::StatusTest::Armijo::check_status(
    NOX::Nln::LineSearch::Controller& ls_controller, const ::NOX::Abstract::Group& grp,
    ::NOX::StatusTest::CheckType checkType)
{
  // setup for the current line search loop
  if (ls_controller.get_num_iterations() == 0)
  {
    // If the search direction is no descent direction,
    // this function detects it and returns the corresponding
    // status.
    if (setup(ls_controller, grp))
      status_ = status_unevaluated;
    else
      status_ = status_no_descent_direction;
  }
  else if (checkType == ::NOX::StatusTest::None)
  {
    fcurr_ = 0.0;
    status_ = status_unevaluated;
  }
  // If the setup call detected a non-descent direction,
  // we will not check the Armijo rule, since the check will
  // fail anyway.
  else if (status_ != status_no_descent_direction)
  {
    const ::NOX::MeritFunction::Generic& mrtFct = ls_controller.get_merit_function();

    fcurr_ = mrtFct.computef(grp);

    step_ = ls_controller.get_step_length();

    // check the Armijo rule
    status_ = (fcurr_ < fref_ + c_1_ * step_ * slope_) ? status_converged : status_step_too_long;
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::StatusType NOX::Nln::Inner::StatusTest::Armijo::get_status() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::Nln::Inner::StatusTest::Armijo::print(std::ostream& stream, int indent) const
{
  std::string indent_string;
  indent_string.assign(indent, ' ');

  stream << indent_string;
  stream << status_;
  if (isMonotone_)
    stream << "Monotone";
  else
    stream << "Non-monotone";
  stream << " ";
  stream << "Armijo-Rule: ";
  stream << ::NOX::Utils::sciformat(fcurr_, 3) << " < "
         << ::NOX::Utils::sciformat(fref_ + c_1_ * step_ * slope_, 3) << "\n";

  stream << indent_string;
  stream << std::setw(13) << " ";
  stream << "(step = " << ::NOX::Utils::sciformat(step_, 3);
  stream << ", slope = " << ::NOX::Utils::sciformat(slope_, 3);
  if (not isMonotone_) stream << ", history = " << maxHistSize_;
  stream << ")\n";

  return stream;
}

FOUR_C_NAMESPACE_CLOSE
