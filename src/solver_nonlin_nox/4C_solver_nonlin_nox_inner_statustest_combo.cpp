/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_inner_statustest_combo.hpp"

#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

// tolerated unconverged states
const std::set<NOX::Nln::Inner::StatusTest::StatusType>
    NOX::Nln::Inner::StatusTest::Combo::unconverged_ = {
        status_step_too_long, status_step_too_short};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::Combo::Combo(
    ::NOX::StatusTest::Combo::ComboType t, const ::NOX::Utils* u)
    : type_(t)
{
  if (u != nullptr) utils_ = *u;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::Combo::Combo(
    ::NOX::StatusTest::Combo::ComboType t, const Teuchos::RCP<Generic>& a, const ::NOX::Utils* u)
    : type_(t)
{
  if (u != nullptr) utils_ = *u;
  // fill ghost vector
  tests_.push_back(a);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::Combo::Combo(::NOX::StatusTest::Combo::ComboType t,
    const Teuchos::RCP<Generic>& a, const Teuchos::RCP<Generic>& b, const ::NOX::Utils* u)
    : type_(t)
{
  if (u != nullptr) utils_ = *u;

  // fill ghost vector
  tests_.push_back(a);
  // Be careful, because the test b was already added to
  // the base class tests vector during the construction call!
  this->addStatusTest(b);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::StatusType NOX::Nln::Inner::StatusTest::Combo::CheckStatus(
    const Interface::Required& interface, const ::NOX::Solver::Generic& solver,
    const ::NOX::Abstract::Group& grp, ::NOX::StatusTest::CheckType checkType)
{
  if (type_ == ::NOX::StatusTest::Combo::OR)
    or_op(interface, solver, grp, checkType);
  else
    and_op(interface, solver, grp, checkType);

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::StatusType NOX::Nln::Inner::StatusTest::Combo::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Inner::StatusTest::Combo::or_op(const Interface::Required& interface,
    const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
    ::NOX::StatusTest::CheckType checkType)
{
  if (checkType == ::NOX::StatusTest::None)
    status_ = status_unevaluated;
  else
    status_ = status_step_too_long;

  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (const auto& test : tests_)
  {
    NOX::Nln::Inner::StatusTest::StatusType s =
        test->CheckStatus(interface, solver, grp, checkType);

    if (unconverged_.find(status_) != unconverged_.end() and
        unconverged_.find(s) == unconverged_.end())
    {
      status_ = s;

      // Turn off checking for the remaining tests
      if (checkType == ::NOX::StatusTest::Minimal) checkType = ::NOX::StatusTest::None;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Inner::StatusTest::Combo::and_op(const Interface::Required& interface,
    const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
    ::NOX::StatusTest::CheckType checkType)
{
  if (checkType == ::NOX::StatusTest::None)
    status_ = status_unevaluated;
  else
    status_ = status_step_too_long;

  bool isUnconverged = false;

  for (const auto& test : tests_)
  {
    NOX::Nln::Inner::StatusTest::StatusType s =
        test->CheckStatus(interface, solver, grp, checkType);

    // If any of the tests are unconverged, then the AND test is
    // unconverged.
    if (unconverged_.find(s) != unconverged_.end())
    {
      isUnconverged = true;
      status_ = s;

      // Turn off checking for the remaining tests
      if (checkType == ::NOX::StatusTest::Minimal) checkType = ::NOX::StatusTest::None;
    }

    // If this is the first test and it's converged/failed, copy its
    // status to the combo status.
    if ((!isUnconverged) and unconverged_.find(status_) != unconverged_.end())
    {
      status_ = s;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::Combo& NOX::Nln::Inner::StatusTest::Combo::addStatusTest(
    const Teuchos::RCP<Generic>& a)
{
  if (is_safe(*a))
  {
    tests_.push_back(a);
  }
  else
  {
    const int indent = 2;
    utils_.err() << "\n*** WARNING! ***\n";
    utils_.err() << "This combo test currently consists of the following:\n";
    this->print(utils_.err(), indent);
    utils_.err() << "Unable to add the following test:\n";
    a->print(utils_.err(), indent);
    utils_.err() << "\n";
  }

  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Inner::StatusTest::Combo::is_safe(Generic& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (&a == this) return false;

  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (auto& test : tests_)
  {
    NOX::Nln::Inner::StatusTest::Combo* ptr =
        dynamic_cast<NOX::Nln::Inner::StatusTest::Combo*>(test.get());
    if (ptr != nullptr)
      if (!ptr->is_safe(a)) return false;
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>>&
NOX::Nln::Inner::StatusTest::Combo::GetTestVector() const
{
  return tests_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::Nln::Inner::StatusTest::Combo::print(std::ostream& stream, int indent) const
{
  stream << std::string(indent, ' ');
  stream << status_;
  stream << ((type_ == ::NOX::StatusTest::Combo::OR) ? "OR" : "AND");
  stream << " Combination";
  stream << " -> " << std::endl;

  for (const auto& test : tests_) test->print(stream, indent + 2);

  return stream;
}

FOUR_C_NAMESPACE_CLOSE
