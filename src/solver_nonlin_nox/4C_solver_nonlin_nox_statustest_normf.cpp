/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN implementation of a NormF status test. This
       test can be used to check the residual (right-hand-side)
       for convergence.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_statustest_normf.hpp"  // class definition

#include "4C_solver_nonlin_nox_group.hpp"

#include <NOX_Solver_Generic.H>
#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::StatusTest::NormF::NormF(const std::vector<NOX::Nln::StatusTest::QuantityType>& checkList,
    const std::vector<::NOX::StatusTest::NormF::ToleranceType>& toltype,
    const std::vector<double>& tolerance,
    const std::vector<::NOX::Abstract::Vector::NormType>& ntype,
    const std::vector<::NOX::StatusTest::NormF::ScaleType>& stype, const ::NOX::Utils* u)
    : nChecks_(checkList.size()),
      checkList_(checkList),
      gStatus_(::NOX::StatusTest::Unevaluated),
      status_(std::vector<::NOX::StatusTest::StatusType>(nChecks_, gStatus_)),
      normType_(ntype),
      scaleType_(stype),
      toleranceType_(toltype),
      specifiedTolerance_(tolerance),
      initialTolerance_(Teuchos::null),
      trueTolerance_(tolerance),
      normF_(Teuchos::null)
{
  if (u != nullptr) utils_ = *u;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::StatusTest::NormF::relative_setup(
    Teuchos::RCP<const ::NOX::Abstract::Group>& initialGuess)
{
  initialTolerance_ = compute_norm(initialGuess);

  if (initialTolerance_.is_null())
  {
    utils_.err() << "NOX::Nln::StatusTest::NormF::NormF - F was not computed for the"
                    " given nox_abstract_group!"
                 << std::endl;
    throw "NOX Error";
  }

  for (std::size_t i = 0; i < nChecks_; ++i)
  {
    if (toleranceType_[i] == ::NOX::StatusTest::NormF::Relative)
      trueTolerance_[i] = specifiedTolerance_[i] * (*initialTolerance_)[i];
    else if (toleranceType_[i] == ::NOX::StatusTest::NormF::Absolute)
      trueTolerance_[i] = specifiedTolerance_[i];
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const std::vector<double>> NOX::Nln::StatusTest::NormF::compute_norm(
    Teuchos::RCP<const ::NOX::Abstract::Group>& grp)
{
  if (!grp->isF()) return Teuchos::null;

  // cast the given group
  Teuchos::RCP<const NOX::Nln::Group> nlnGrp =
      Teuchos::rcp_dynamic_cast<const NOX::Nln::Group>(grp);
  if (nlnGrp.is_null())
  {
    utils_.err() << "::NOX::StatusTest::NormF::ComputeNorm - NOX::Nln::Group cast failed!"
                 << std::endl;
    throw "NOX Error";
  }

  Teuchos::RCP<const std::vector<double>> norms =
      nlnGrp->get_rhs_norms(normType_, checkList_, Teuchos::rcp(&scaleType_, false));

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::Nln::StatusTest::NormF::checkStatus(
    const ::NOX::Solver::Generic& problem, ::NOX::StatusTest::CheckType checkType)
{
  if (problem.getNumIterations() == 0)
  {
    Teuchos::RCP<const ::NOX::Abstract::Group> grp =
        Teuchos::rcpFromRef(problem.getSolutionGroup());
    relative_setup(grp);
  }

  if (checkType == ::NOX::StatusTest::None)
  {
    normF_ = Teuchos::rcp(new std::vector<double>(nChecks_, 0.0));
    status_.assign(status_.size(), ::NOX::StatusTest::Unevaluated);
    gStatus_ = ::NOX::StatusTest::Unevaluated;
  }
  else
  {
    Teuchos::RCP<const ::NOX::Abstract::Group> grp =
        Teuchos::rcpFromRef(problem.getSolutionGroup());
    normF_ = compute_norm(grp);

    gStatus_ = ::NOX::StatusTest::Converged;
    for (std::size_t i = 0; i < normF_->size(); ++i)
    {
      status_[i] = ((not normF_.is_null()) && ((*normF_)[i] < trueTolerance_[i]))
                       ? ::NOX::StatusTest::Converged
                       : ::NOX::StatusTest::Unconverged;
      if (status_[i] == ::NOX::StatusTest::Unconverged) gStatus_ = ::NOX::StatusTest::Unconverged;
    }
  }

  return gStatus_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::StatusTest::NormF::get_norm_f(
    const NOX::Nln::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < nChecks_; ++i)
    if (checkList_[i] == qType) return (*normF_)[i];

  // if we cannot find the right quantity in the class list we will return -1.0
  return (-1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::StatusTest::NormF::get_true_tolerance(
    const NOX::Nln::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < nChecks_; ++i)
    if (checkList_[i] == qType) return trueTolerance_[i];

  // if we cannot find the right quantity in the class list we will return -1.0
  return (-1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::StatusTest::NormF::get_specified_tolerance(
    const NOX::Nln::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < nChecks_; ++i)
    if (checkList_[i] == qType) return specifiedTolerance_[i];

  // if we cannot find the right quantity in the class list we will return -1.0
  return (-1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::StatusTest::NormF::GetInitialTolerance(
    const NOX::Nln::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < nChecks_; ++i)
    if (checkList_[i] == qType) return (*initialTolerance_)[i];

  // if we cannot find the right quantity in the class list we will return -1.0
  return (-1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::Nln::StatusTest::NormF::get_norm_type(
    const NOX::Nln::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < nChecks_; ++i)
    if (checkList_[i] == qType) return normType_[i];

  // if we cannot find the right quantity in the class list we will return -100
  return (-100);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::StatusTest::NormF::is_quantity(const NOX::Nln::StatusTest::QuantityType& qType) const
{
  for (std::size_t i = 0; i < nChecks_; ++i)
    if (checkList_[i] == qType) return true;

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::Nln::StatusTest::NormF::getStatus() const { return gStatus_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::Nln::StatusTest::NormF::print(std::ostream& stream, int indent) const
{
  std::string indent_string;
  indent_string.assign(indent, ' ');

  for (std::size_t i = 0; i < nChecks_; ++i)
  {
    stream << indent_string;
    stream << status_[i];
    stream << QuantityType2String(checkList_[i]) << "-";
    stream << "F-Norm = " << ::NOX::Utils::sciformat((*normF_)[i], 3);
    stream << " < " << ::NOX::Utils::sciformat(trueTolerance_[i], 3);
    stream << std::endl;

    stream << indent_string;
    stream << std::setw(13) << " ";
    stream << "(";

    if (scaleType_[i] == ::NOX::StatusTest::NormF::Scaled)
      stream << "Scaled";
    else
      stream << "Unscaled";

    stream << " ";

    if (normType_[i] == ::NOX::Abstract::Vector::TwoNorm)
      stream << "2-Norm";
    else if (normType_[i] == ::NOX::Abstract::Vector::OneNorm)
      stream << "1-Norm";
    else if (normType_[i] == ::NOX::Abstract::Vector::MaxNorm)
      stream << "Max-Norm";

    stream << ", ";

    if (toleranceType_[i] == ::NOX::StatusTest::NormF::Absolute)
      stream << "Absolute";
    else
      stream << "Relative";

    stream << ")";

    stream << std::endl;
  }
  return stream;
}

FOUR_C_NAMESPACE_CLOSE
