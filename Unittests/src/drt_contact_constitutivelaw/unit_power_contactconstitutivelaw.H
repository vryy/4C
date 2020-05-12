/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the power constact constitutivelaw

\level 2

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/
#ifndef UNIT_POWER_COCONSTLAW_H
#define UNIT_POWER_COCONSTLAW_H

#include <cxxtest/TestSuite.h>

#include "../../../src/drt_contact_constitutivelaw/power_contactconstitutivelaw.H"

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class PowerConstitutiveLaw_TestSuite;
  }
}  // namespace CONTACT

// class implementation
class CONTACT::CONSTITUTIVELAW::PowerConstitutiveLaw_TestSuite : public CxxTest::TestSuite
{
 public:
  //! set up unit testing
  void setUp()
  {
    // initialize container for material parameters
    const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Container(
            1, INPAR::CONTACT::ConstitutiveLawType::colaw_power, "Power Constitutivelaw"));

    // add parameters to container
    container->Add("A", 3.0);
    container->Add("B", 3.0);
    container->Add("Offset", 0.5);

    const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> powercoconstlaw =
        CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(container);
    coconstlaw_ = powercoconstlaw;

    return;
  }

  //! tear down unit testing
  void tearDown()
  {
    coconstlaw_ = Teuchos::null;
    return;
  }

  //! test member function Evaluate
  void test_Evaluate()
  {
    // gap < 0
    TS_ASSERT_THROWS_ANYTHING(coconstlaw_->Evaluate(1.0));
    // 0< gap < offset
    TS_ASSERT_THROWS_ANYTHING(coconstlaw_->Evaluate(-0.25));
    // offset < gap
    TS_ASSERT_DELTA(coconstlaw_->Evaluate(-0.75), -0.046875, 1.e-15);

    return;
  }

  //! test member function Evaluate
  void test_EvaluateDeriv()
  {
    TS_ASSERT_DELTA(coconstlaw_->EvaluateDeriv(-0.75), 0.5625, 1.e-15);
    TS_ASSERT_THROWS_ANYTHING(coconstlaw_->EvaluateDeriv(-0.25));
    return;
  }

 private:
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;
};
#endif
