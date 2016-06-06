/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_normwrms.cpp

\brief %NOX::NLN Weighted root mean square test of the solution
       increment. A detailed description can be found in the NOX
       documentation.

\maintainer Michael Hiermeier

\date Aug 4, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_statustest_normwrms.H"
#include "nox_nln_group.H"

#include <NOX_Solver_LineSearchBased.H>
#include <NOX_Utils.H>

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::NormWRMS::NormWRMS(
    const std::vector<NOX::NLN::StatusTest::QuantityType>& checkList,
    const std::vector<double>& rtol,
    const std::vector<double>& atol,
    const std::vector<double>& BDFMultiplier,
    const std::vector<double>& tolerance,
    const double& alpha,
    const double& beta,
    const std::vector<bool>& disable_implicit_weighting) :
    normWRMS_(Teuchos::null),
    nChecks_(checkList.size()),
    checkList_(checkList),
    rtol_(rtol),
    atol_(atol),
    factor_(BDFMultiplier),
    tol_(tolerance),
    alpha_(alpha),
    computedStepSize_(1.0),
    beta_(beta),
    achievedTol_(0.0),
    gStatus_(NOX::StatusTest::Unconverged),
    status_(std::vector<NOX::StatusTest::StatusType>(nChecks_,gStatus_)),
    printCriteria2Info_(false),
    printCriteria3Info_(false),
    disable_implicit_weighting_(disable_implicit_weighting)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::StatusTest::NormWRMS::checkStatus(
    const NOX::Solver::Generic& problem,
    NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
  {
    gStatus_ = NOX::StatusTest::Unevaluated;
    status_.assign(nChecks_,gStatus_);
    normWRMS_ = Teuchos::rcp(new std::vector<double>(nChecks_,1.0e+12));
    return gStatus_;
  }

  Teuchos::RCP<const Abstract::Group> soln =
      Teuchos::rcpFromRef(problem.getSolutionGroup());
  // all entries of the status_ vector are initialized to a unconverged status
  gStatus_ = NOX::StatusTest::Unconverged;
  status_ = std::vector<NOX::StatusTest::StatusType>(nChecks_,gStatus_);

  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0)
  {
    normWRMS_ = Teuchos::rcp(new std::vector<double>(nChecks_,1.0e+12));
    return gStatus_;
  }

  // ---------------------------------------------------------
  // Begin check for convergence criteria #1 (local check)
  // ---------------------------------------------------------
  // cast the nox_abstract_group to nox_nln_group
  Teuchos::RCP<const NOX::NLN::Group> nlnGrp =
      Teuchos::rcp_dynamic_cast<const NOX::NLN::Group>(soln);

  // all entries of the criteria vector are initialized to Converged status
  std::vector<NOX::StatusTest::StatusType> criteria =
      std::vector<NOX::StatusTest::StatusType>(3,NOX::StatusTest::Converged);

  // get the solution vector of the last step
  const NOX::Abstract::Vector& xOld =
      problem.getPreviousSolutionGroup().getX();

  // get the root mean square from the underlying interface classes
  normWRMS_ = nlnGrp->GetSolutionUpdateRMS(xOld,atol_,rtol_,checkList_,disable_implicit_weighting_);

  // loop over all quantities
  for (std::size_t i=0;i<nChecks_;++i)
  {
    // do the weighting by the given factor
    normWRMS_->at(i) *= factor_.at(i);
    if (normWRMS_->at(i) <= tol_.at(i))
      status_.at(i) = NOX::StatusTest::Converged;
    else
      criteria[0] = NOX::StatusTest::Unconverged;
  }

  // ---------------------------------------------------------
  // Begin check for convergence criteria #2 (global check)
  // ---------------------------------------------------------
  // Determine if the Generic solver is a LineSearchBased solver
  // If it is not then return a "Converged" status
  const NOX::Solver::Generic* test = NULL;
  test = dynamic_cast<const NOX::Solver::LineSearchBased*>(&problem);
  if (test == NULL)
    criteria[1] = NOX::StatusTest::Converged;
  else
  {
    printCriteria2Info_ = true;
    computedStepSize_ =
      (dynamic_cast<const NOX::Solver::LineSearchBased*>(&problem))->getStepSize();

    if (computedStepSize_ < alpha_)
    {
      status_.assign(nChecks_,NOX::StatusTest::Unconverged);
      criteria[1] = NOX::StatusTest::Unconverged;
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
            const Teuchos::ParameterList& list = p.sublist("Direction").
                sublist("Newton").sublist("Linear Solver").sublist("Output");
            if (Teuchos::isParameterType<double>(list, "Achieved Tolerance"))
              printCriteria3Info_ = true;
          }

  if (printCriteria3Info_)
  {
    achievedTol_ = const_cast<Teuchos::ParameterList&>(problem.getList()).
      sublist("Direction").sublist("Newton").sublist("Linear Solver").
      sublist("Output").get("Achieved Tolerance", -1.0);
    if (achievedTol_ > beta_)
    {
      criteria[2] =  NOX::StatusTest::Unconverged;
      status_.assign(nChecks_,NOX::StatusTest::Unconverged);
    }
  }

  // Determine global status of test
  gStatus_ =   ((criteria[0] == NOX::StatusTest::Converged) and
                (criteria[1] == NOX::StatusTest::Converged) and
                (criteria[2] == NOX::StatusTest::Converged))
                ? NOX::StatusTest::Converged : NOX::StatusTest::Unconverged;

  return gStatus_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::StatusTest::NormWRMS::IsQuantity(
    const NOX::NLN::StatusTest::QuantityType& qType) const
{
  for (std::size_t i=0;i<nChecks_;++i)
    if (checkList_[i]==qType)
      return true;

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::StatusTest::NormWRMS::GetAbsoluteTolerance(
    const NOX::NLN::StatusTest::QuantityType& qType) const
{
  for (std::size_t i=0;i<nChecks_;++i)
    if (checkList_[i]==qType)
      return atol_[i];

  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::StatusTest::NormWRMS::GetRelativeTolerance(
    const NOX::NLN::StatusTest::QuantityType& qType) const
{
  for (std::size_t i=0;i<nChecks_;++i)
    if (checkList_[i]==qType)
      return rtol_[i];

  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::StatusTest::NormWRMS::getStatus() const
{
  return gStatus_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::StatusTest::NormWRMS::print(
    std::ostream& stream,
    int indent) const
{
  std::string indent_string;
  indent_string.assign(indent,' ');

  for (std::size_t i=0;i<nChecks_;++i)
  {
    stream << indent_string;
    stream << status_[i];
    stream << QuantityType2String(checkList_[i]) << "-";
    stream << "WRMS-Norm = " << NOX::Utils::sciformat((*normWRMS_)[i], 3)
           << " < " << tol_[i];
    stream << std::endl;
  }
  if (printCriteria2Info_) {
    stream << indent_string;
    stream << std::setw(13) << " ";
    stream << "(Min Step Size:  " << NOX::Utils::sciformat(computedStepSize_, 3)
           << " >= " << alpha_ << ")";
    stream << std::endl;
  }
  if (printCriteria3Info_) {
    stream << indent_string;
    stream << std::setw(13) << " ";
    stream << "(Max Lin Solv Tol:  " << NOX::Utils::sciformat(achievedTol_, 3)
           << " < " << beta_ << ")";
    stream << std::endl;
  }
  return stream;
}
