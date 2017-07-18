/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_normupdate.cpp

\brief %NOX::NLN implementation of a NormUpdate status test. This
       test can be used to check the solution increment \f$\Delta x\f$
       for convergence.

\maintainer Michael Hiermeier

\date Sep 17, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_statustest_normupdate.H"
#include "nox_nln_group.H"

#include <NOX_Solver_LineSearchBased.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::NormUpdate::NormUpdate(
    const std::vector<NOX::NLN::StatusTest::QuantityType>& checkList,
    const std::vector<NormUpdate::ToleranceType>& toltype,
    const std::vector<double>& tolerance,
    const std::vector<NOX::Abstract::Vector::NormType>& ntype,
    const std::vector<NormUpdate::ScaleType>& stype,
    const double& alpha,
    const double& beta,
    const NOX::Utils* u)
    : nChecks_(checkList.size()),
      computedStepSize_(1.0),
      alpha_(alpha),
      achievedTol_(0.0),
      beta_(beta),
      checkList_(checkList),
      gStatus_(NOX::StatusTest::Unevaluated),
      status_(std::vector<NOX::StatusTest::StatusType>(nChecks_,gStatus_)),
      normType_(ntype),
      scaleType_(stype),
      toleranceType_(toltype),
      specifiedTolerance_(tolerance),
      normRefSol_(Teuchos::null),
      trueTolerance_(tolerance),
      normUpdate_(Teuchos::null),
      printCriteria2Info_(false),
      printCriteria3Info_(false)
{
  if (u!=NULL)
    utils_ = *u;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::StatusTest::NormUpdate::ComputeNorm(
    const NOX::Abstract::Group& grp,
    const NOX::Solver::Generic& problem)
{
  // cast the nox_abstract_group to nox_nln_group
  const NOX::NLN::Group* nlngrp = dynamic_cast<const NOX::NLN::Group*>(&grp);
  if (nlngrp == NULL)
    throwError("ComputeNorm","Dynamic cast to NOX::NLN::Group failed!");

  // get the old solution vector
  const NOX::Abstract::Vector& xOld = problem.getPreviousSolutionGroup().getX();

  // (1) of the increment of the given quantities
  normUpdate_ = nlngrp->GetSolutionUpdateNorms(
      xOld,normType_,checkList_,Teuchos::rcp(&scaleType_,false));
  // (2) of the last accepted Newton step
  normRefSol_ = nlngrp->GetPreviousSolutionNorms(
      xOld,normType_,checkList_,Teuchos::rcp(&scaleType_,false));

  for (std::size_t i=0;i<nChecks_;++i)
  {
    if (toleranceType_[i]==NormUpdate::Relative)
      trueTolerance_[i] = specifiedTolerance_[i] * normRefSol_->at(i);
    else if (toleranceType_[i]==NormUpdate::Absolute)
      trueTolerance_[i] = specifiedTolerance_[i];
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::StatusTest::NormUpdate::checkStatus(
    const NOX::Solver::Generic& problem,
    NOX::StatusTest::CheckType checkType
    )
{
  if (checkType == NOX::StatusTest::None)
  {
    gStatus_ = NOX::StatusTest::Unevaluated;
    status_.assign(nChecks_,gStatus_);
    normUpdate_ = Teuchos::rcp(new std::vector<double>(nChecks_,1.0e+12));
    normRefSol_ = Teuchos::rcp(new std::vector<double>(nChecks_,1.0));
    return gStatus_;
  }

  const NOX::Abstract::Group& soln = problem.getSolutionGroup();
  // all entries of the status_ vector are initialized to a unconverged status
  gStatus_ = NOX::StatusTest::Unconverged;
  status_ = std::vector<NOX::StatusTest::StatusType>(nChecks_,gStatus_);

  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0)
  {
    // set some default values
    normUpdate_ = Teuchos::rcp(new std::vector<double>(nChecks_,1.0e+12));
    normRefSol_ = Teuchos::rcp(new std::vector<double>(nChecks_,1.0));
    return gStatus_;
  }

  // ---------------------------------------------------------
  // Begin check for convergence criteria #1 (local check)
  // ---------------------------------------------------------
  // all entries of the criteria vector are initialized to Converged status
  std::vector<NOX::StatusTest::StatusType> criteria =
      std::vector<NOX::StatusTest::StatusType>(3,NOX::StatusTest::Converged);

  // get the specified norms from the underlying interface classes.
  // update of the truetolerance variable.
  ComputeNorm(soln,problem);

  // loop over all quantities
  for (std::size_t i=0;i<nChecks_;++i)
  {
    if (normUpdate_->at(i) <= trueTolerance_.at(i))
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
  // First time through, make sure the output parameter list exists.
  // Since the list is const, a sublist call to a non-existent sublist
  // throws an error.  Therefore we have to check the existence of each
  // sublist before we call it.
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
int NOX::NLN::StatusTest::NormUpdate::GetNormType(
    const NOX::NLN::StatusTest::QuantityType& qType) const
{
  for (std::size_t i=0;i<nChecks_;++i)
    if (checkList_[i]==qType)
      return normType_[i];

  // if we cannot find the right quantity in the class list we will return -100
  return (-100);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::StatusTest::NormUpdate::IsQuantity(
    const NOX::NLN::StatusTest::QuantityType& qType) const
{
  for (std::size_t i=0;i<nChecks_;++i)
    if (checkList_[i]==qType)
      return true;

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::StatusTest::NormUpdate::getStatus() const
{
  return gStatus_;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::StatusTest::NormUpdate::print(
    std::ostream& stream,
    int indent) const
{
  std::string indent_string;
  indent_string.assign(indent,' ');

  for (std::size_t i=0;i<nChecks_;++i)
  {
    if (i>0) stream << "\n";
    stream << indent_string;
    stream << status_[i];
    stream << QuantityType2String(checkList_[i]) << "-";
    stream << "Update-Norm = " << NOX::Utils::sciformat(normUpdate_->at(i), 3)
           << " < " << trueTolerance_[i];
    stream << std::endl;

    stream << indent_string;
    stream << std::setw(13) << " ";
    stream << "(";

    if (scaleType_[i] == NormUpdate::Scaled)
      stream << "Scaled";
    else
      stream << "Unscaled";

    stream << " ";

    if (normType_[i] == NOX::Abstract::Vector::TwoNorm)
      stream << "2-Norm";
    else if (normType_[i] == NOX::Abstract::Vector::OneNorm)
      stream << "1-Norm";
    else if (normType_[i] == NOX::Abstract::Vector::MaxNorm)
      stream << "Max-Norm";

    stream << ", ";

    if (toleranceType_[i] == NormUpdate::Absolute)
      stream << "Absolute";
    else
      stream << "Relative";

    stream << ")";
  }

  if (printCriteria2Info_) {
    stream << std::endl << indent_string;
    stream << std::setw(13) << " ";
    stream << "(Min Step Size:  " << NOX::Utils::sciformat(computedStepSize_, 3)
           << " >= " << alpha_ << ")";
  }
  if (printCriteria3Info_) {
    stream << std::endl << indent_string;
    stream << std::setw(13) << " ";
    stream << "(Max Lin Solv Tol:  " << NOX::Utils::sciformat(achievedTol_, 3)
           << " < " << beta_ << ")";
  }
  stream << std::endl;

  return stream;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::StatusTest::NormUpdate::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error)) {
    utils_.out() << "ERROR - NOX::NLN::StatusTest::NormUpdate::" << functionName
     << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}
