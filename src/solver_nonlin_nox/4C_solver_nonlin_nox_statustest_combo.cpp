/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_statustest_combo.hpp"

#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::StatusTest::Combo::Combo(ComboType t, const ::NOX::Utils* u)
    : ::NOX::StatusTest::Combo(t, u)
{
  if (u != nullptr) utils_ = *u;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::StatusTest::Combo::Combo(
    ComboType t, const Teuchos::RCP<Generic>& a, const ::NOX::Utils* u)
    : ::NOX::StatusTest::Combo(t, a, u)
{
  if (u != nullptr) utils_ = *u;
  // fill ghost vector
  tests_.push_back(a);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::StatusTest::Combo::Combo(ComboType t, const Teuchos::RCP<Generic>& a,
    const Teuchos::RCP<Generic>& b, const ::NOX::Utils* u)
    : ::NOX::StatusTest::Combo(t, a, b, u)
{
  if (u != nullptr) utils_ = *u;

  // fill ghost vector
  tests_.push_back(a);
  // Be careful, because the test b was already added to
  // the base class tests vector during the construction call!
  this->add_status_test(b, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::StatusTest::Combo& NOX::Nln::StatusTest::Combo::addStatusTest(
    const Teuchos::RCP<::NOX::StatusTest::Generic>& a)
{
  return add_status_test(a, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::StatusTest::Combo& NOX::Nln::StatusTest::Combo::add_status_test(
    const Teuchos::RCP<::NOX::StatusTest::Generic>& a, const bool& init)
{
  if (isSafe(*(a.get())))
  {
    tests_.push_back(a);
    // add the test to the test-vector of the base class
    if (not init) ::NOX::StatusTest::Combo::addStatusTest(a);
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
bool NOX::Nln::StatusTest::Combo::isSafe(::NOX::StatusTest::Generic& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (&a == this) return false;

  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (std::vector<Teuchos::RCP<::NOX::StatusTest::Generic>>::iterator i = tests_.begin();
       i != tests_.end(); ++i)
  {
    NOX::Nln::StatusTest::Combo* ptr = dynamic_cast<NOX::Nln::StatusTest::Combo*>(i->get());
    if (ptr != nullptr)
      if (!ptr->isSafe(a)) return false;
  }

  // call base version
  return ::NOX::StatusTest::Combo::isSafe(a);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<Teuchos::RCP<::NOX::StatusTest::Generic>>&
NOX::Nln::StatusTest::Combo::get_test_vector() const
{
  return tests_;
}

FOUR_C_NAMESPACE_CLOSE
