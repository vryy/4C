// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_inner_statustest_upperbound.hpp"

#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_linesearch_controller.hpp"

#include <NOX_Abstract_Vector.H>
#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::UpperBound::UpperBound(const double& upperboundval,
    const ::NOX::Abstract::Vector::NormType normtype,
    const NOX::Nln::StatusTest::QuantityType qtype)
    : status_(status_unevaluated),
      normtype_(normtype),
      qtype_(qtype),
      upperboundval_(upperboundval),
      reduction_fac_(1.0),
      stepmaxval_(0.0)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Inner::StatusTest::UpperBound::get_search_direction_length(
    const NOX::Nln::LineSearch::Controller& ls_controller, const ::NOX::Abstract::Group& grp) const
{
  const NOX::Nln::Group& nln_grp = dynamic_cast<const NOX::Nln::Group&>(grp);

  return nln_grp.get_trial_update_norm(ls_controller.get_search_direction(), normtype_, qtype_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::StatusType NOX::Nln::Inner::StatusTest::UpperBound::check_status(
    NOX::Nln::LineSearch::Controller& ls_controller, const ::NOX::Abstract::Group& grp,
    ::NOX::StatusTest::CheckType checkType)
{
  /* we reduce the step length according to the upper bound criterion in the first
   * line search (i.e. inner) iteration and do nothing in all following iterations */
  if (ls_controller.get_num_iterations() == 0)
  {
    const double dir_length = get_search_direction_length(ls_controller, grp);
    double steplength = ls_controller.get_step_length();

    // compute specified norm
    stepmaxval_ = steplength * dir_length;

    // check the value for specified upper bound
    status_ = (stepmaxval_ < upperboundval_) ? status_converged : status_step_too_long;

    /* in opposite to classical line search, there is no need to find the new step
     * length iteratively.
     * instead, we can immediately reduce the step length appropriately and avoid
     * several unnecessary evaluations of rhs (computeF calls). */
    if (status_ == status_step_too_long)
    {
      /* the following is equivalent to dividing the step successively by two until
       * criterion is met. note: upperboundval_!=0 is checked in
       * NOX::Nln::Inner::StatusTest::Factory::build_upper_bound_test */
      reduction_fac_ =
          std::pow(0.5, std::ceil(std::log(stepmaxval_ / upperboundval_) / std::log(2)));

      steplength *= reduction_fac_;
      ls_controller.set_step_length(steplength);

      // adapt the stepmaxval_ variable accordingly to get correct output from print()
      stepmaxval_ *= reduction_fac_;

      status_ = status_converged;
    }
    // nothing to do. reset the the reduction factor to its default value
    else
      reduction_fac_ = 1.0;
  }
  else if (checkType == ::NOX::StatusTest::None)
  {
    /*!
      The test can (and should, if possible) be skipped if
      checkType is NOX::StatusType::None. If the test is skipped, then
      the status should be set to ::NOX::StatusTest::Unevaluated.
    */
    status_ = status_unevaluated;
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::StatusType NOX::Nln::Inner::StatusTest::UpperBound::get_status() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::Nln::Inner::StatusTest::UpperBound::print(std::ostream& stream, int indent) const
{
  std::string indent_string;
  indent_string.assign(indent, ' ');

  stream << indent_string;
  stream << status_;
  stream << " ";
  stream << "Upper Bound for Newton step size: ";
  stream << ::NOX::Utils::sciformat(stepmaxval_, 3) << " < "
         << ::NOX::Utils::sciformat(upperboundval_, 3) << "\n";
  stream << indent_string << ::NOX::Utils::fill(13, ' ')
         << " (reduction factor = " << ::NOX::Utils::sciformat(reduction_fac_, 3) << ")\n";

  return stream;
}

FOUR_C_NAMESPACE_CLOSE
