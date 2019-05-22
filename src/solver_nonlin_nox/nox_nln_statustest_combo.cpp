/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_combo.cpp

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_statustest_combo.H"

#include <NOX_Utils.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::Combo::Combo(ComboType t, const NOX::Utils* u) : NOX::StatusTest::Combo(t, u)
{
  if (u != NULL) utils_ = *u;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::Combo::Combo(ComboType t, const Teuchos::RCP<Generic>& a, const NOX::Utils* u)
    : NOX::StatusTest::Combo(t, a, u)
{
  if (u != NULL) utils_ = *u;
  // fill ghost vector
  tests_.push_back(a);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::Combo::Combo(ComboType t, const Teuchos::RCP<Generic>& a,
    const Teuchos::RCP<Generic>& b, const NOX::Utils* u)
    : NOX::StatusTest::Combo(t, a, b, u)
{
  if (u != NULL) utils_ = *u;

  // fill ghost vector
  tests_.push_back(a);
  // Be careful, because the test b was already added to
  // the base class tests vector during the construction call!
  this->addStatusTest(b, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::Combo& NOX::NLN::StatusTest::Combo::addStatusTest(
    const Teuchos::RCP<NOX::StatusTest::Generic>& a)
{
  return addStatusTest(a, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::Combo& NOX::NLN::StatusTest::Combo::addStatusTest(
    const Teuchos::RCP<NOX::StatusTest::Generic>& a, const bool& init)
{
  if (isSafe(*(a.get())))
  {
    tests_.push_back(a);
    // add the test to the test-vector of the base class
    if (not init) NOX::StatusTest::Combo::addStatusTest(a);
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
bool NOX::NLN::StatusTest::Combo::isSafe(NOX::StatusTest::Generic& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (&a == this) return false;

  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (std::vector<Teuchos::RCP<NOX::StatusTest::Generic>>::iterator i = tests_.begin();
       i != tests_.end(); ++i)
  {
    NOX::NLN::StatusTest::Combo* ptr = dynamic_cast<NOX::NLN::StatusTest::Combo*>(i->get());
    if (ptr != NULL)
      if (!ptr->isSafe(a)) return false;
  }

  // call base version
  return NOX::StatusTest::Combo::isSafe(a);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<Teuchos::RCP<NOX::StatusTest::Generic>>&
NOX::NLN::StatusTest::Combo::GetTestVector() const
{
  return tests_;
}
