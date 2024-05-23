/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN Weighted root mean square test of the solution
       increment. A detailed description can be found in the NOX
       documentation.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_statustest_normwrms.hpp"

#include "4C_solver_nonlin_nox_group.hpp"

#include <NOX_Solver_LineSearchBased.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::NormWRMS::NormWRMS(
    const std::vector<NOX::NLN::StatusTest::QuantityType>& checkList,
    const std::vector<double>& rtol, const std::vector<double>& atol,
    const std::vector<double>& BDFMultiplier, const std::vector<double>& tolerance,
    const double& alpha, const double& beta, const std::vector<bool>& disable_implicit_weighting)
    : norm_wrms_(Teuchos::null),
      n_checks_(checkList.size()),
      check_list_(checkList),
      rtol_(rtol),
      atol_(atol),
      factor_(BDFMultiplier),
      tol_(tolerance),
      alpha_(alpha),
      computed_step_size_(1.0),
      beta_(beta),
      achieved_tol_(0.0),
      g_status_(::NOX::StatusTest::Unconverged),
      status_(std::vector<::NOX::StatusTest::StatusType>(n_checks_, g_status_)),
      print_criteria2_info_(false),
      print_criteria3_info_(false),
      disable_implicit_weighting_(disable_implicit_weighting)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::NLN::StatusTest::NormWRMS::checkStatus(
    const ::NOX::Solver::Generic& problem, ::NOX::StatusTest::CheckType checkType)
{
  if (checkType == ::NOX::StatusTest::None)
  {
    g_status_ = ::NOX::StatusTest::Unevaluated;
    status_.assign(n_checks_, g_status_);
    norm_wrms_ = Teuchos::rcp(new std::vector<double>(n_checks_, 1.0e+12));
    return g_status_;
  }

  Teuchos::RCP<const ::NOX::Abstract::Group> soln = Teuchos::rcpFromRef(problem.getSolutionGroup());
  // all entries of the status_ vector are initialized to a unconverged status
  g_status_ = ::NOX::StatusTest::Unconverged;
  status_ = std::vector<::NOX::StatusTest::StatusType>(n_checks_, g_status_);

  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0)
  {
    norm_wrms_ = Teuchos::rcp(new std::vector<double>(n_checks_, 1.0e+12));
    return g_status_;
  }

  // ---------------------------------------------------------
  // Begin check for convergence criteria #1 (local check)
  // ---------------------------------------------------------
  // cast the nox_abstract_group to nox_nln_group
  Teuchos::RCP<const NOX::NLN::Group> nlnGrp =
      Teuchos::rcp_dynamic_cast<const NOX::NLN::Group>(soln);

  // all entries of the criteria vector are initialized to Converged status
  std::vector<::NOX::StatusTest::StatusType> criteria =
      std::vector<::NOX::StatusTest::StatusType>(3, ::NOX::StatusTest::Converged);

  // get the solution vector of the last step
  const ::NOX::Abstract::Vector& xOld = problem.getPreviousSolutionGroup().getX();

  // get the root mean square from the underlying interface classes
  norm_wrms_ =
      nlnGrp->get_solution_update_rms(xOld, atol_, rtol_, check_list_, disable_implicit_weighting_);

  // loop over all quantities
  for (std::size_t i = 0; i < n_checks_; ++i)
  {
    // do the weighting by the given factor
    norm_wrms_->at(i) *= factor_.at(i);
    if (norm_wrms_->at(i) <= tol_.at(i))
      status_.at(i) = ::NOX::StatusTest::Converged;
    else
      criteria[0] = ::NOX::StatusTest::Unconverged;
  }

  // ---------------------------------------------------------
  // Begin check for convergence criteria #2 (global check)
  // ---------------------------------------------------------
  // Determine if the Generic solver is a LineSearchBased solver
  // If it is not then return a "Converged" status
  const ::NOX::Solver::Generic* test = nullptr;
  test = dynamic_cast<const ::NOX::Solver::LineSearchBased*>(&problem);
  if (test == nullptr)
    criteria[1] = ::NOX::StatusTest::Converged;
  else
  {
    print_criteria2_info_ = true;
    computed_step_size_ =
        (dynamic_cast<const ::NOX::Solver::LineSearchBased*>(&problem))->getStepSize();

    if (computed_step_size_ < alpha_)
    {
      status_.assign(n_checks_, ::NOX::StatusTest::Unconverged);
      criteria[1] = ::NOX::StatusTest::Unconverged;
    }
  }

  // ---------------------------------------------------------
  // Begin check for convergence criteria #3 (global check)
  // ---------------------------------------------------------
  /* First time through, make sure the output parameter list exists.
   * Since the list is const, a sublist call to a non-existent sublist
   * throws an error.  Therefore we have to check the existence of each
   * sublist before we call it. */
  const Teuchos::ParameterList& p = problem.getList();
  if (niters == 1)
    if (p.isSublist("Direction"))
      if (p.sublist("Direction").isSublist("Newton"))
        if (p.sublist("Direction").sublist("Newton").isSublist("Linear Solver"))
          if (p.sublist("Direction").sublist("Newton").sublist("Linear Solver").isSublist("Output"))
          {
            const Teuchos::ParameterList& list =
                p.sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output");
            if (Teuchos::isParameterType<double>(list, "Achieved Tolerance"))
              print_criteria3_info_ = true;
          }

  if (print_criteria3_info_)
  {
    achieved_tol_ = const_cast<Teuchos::ParameterList&>(problem.getList())
                        .sublist("Direction")
                        .sublist("Newton")
                        .sublist("Linear Solver")
                        .sublist("Output")
                        .get("Achieved Tolerance", -1.0);
    if (achieved_tol_ > beta_)
    {
      criteria[2] = ::NOX::StatusTest::Unconverged;
      status_.assign(n_checks_, ::NOX::StatusTest::Unconverged);
    }
  }

  // Determine global status of test
  g_status_ = ((criteria[0] == ::NOX::StatusTest::Converged) and
                  (criteria[1] == ::NOX::StatusTest::Converged) and
                  (criteria[2] == ::NOX::StatusTest::Converged))
                  ? ::NOX::StatusTest::Converged
                  : ::NOX::StatusTest::Unconverged;

  return g_status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::StatusTest::NormWRMS::IsQuantity(
    const NOX::NLN::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < n_checks_; ++i)
    if (check_list_[i] == qType) return true;

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::StatusTest::NormWRMS::get_absolute_tolerance(
    const NOX::NLN::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < n_checks_; ++i)
    if (check_list_[i] == qType) return atol_[i];

  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::StatusTest::NormWRMS::get_relative_tolerance(
    const NOX::NLN::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < n_checks_; ++i)
    if (check_list_[i] == qType) return rtol_[i];

  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::NLN::StatusTest::NormWRMS::getStatus() const
{
  return g_status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::StatusTest::NormWRMS::print(std::ostream& stream, int indent) const
{
  std::string indent_string;
  indent_string.assign(indent, ' ');

  for (std::size_t i = 0; i < n_checks_; ++i)
  {
    stream << indent_string;
    stream << status_[i];
    stream << QuantityType2String(check_list_[i]) << "-";
    stream << "WRMS-Norm = " << ::NOX::Utils::sciformat((*norm_wrms_)[i], 3) << " < "
           << ::NOX::Utils::sciformat(tol_[i], 3);
    stream << std::endl;
  }
  if (print_criteria2_info_)
  {
    stream << indent_string;
    stream << std::setw(13) << " ";
    stream << "(Min Step Size:  " << ::NOX::Utils::sciformat(computed_step_size_, 3)
           << " >= " << ::NOX::Utils::sciformat(alpha_, 3) << ")";
    stream << std::endl;
  }
  if (print_criteria3_info_)
  {
    stream << indent_string;
    stream << std::setw(13) << " ";
    stream << "(Max Lin Solv Tol:  " << ::NOX::Utils::sciformat(achieved_tol_, 3) << " < "
           << ::NOX::Utils::sciformat(beta_, 3) << ")";
    stream << std::endl;
  }
  return stream;
}

FOUR_C_NAMESPACE_CLOSE
