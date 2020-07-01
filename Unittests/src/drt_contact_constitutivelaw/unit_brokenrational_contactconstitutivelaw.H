/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the brokenrational contact constitutivelaw

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef UNIT_BROKENRATIONAL_COCONSTLAW_H
#define UNIT_BROKENRATIONAL_COCONSTLAW_H

#include <cxxtest/TestSuite.h>

#include "src/drt_contact_constitutivelaw/brokenrational_contactconstitutivelaw.H"

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class BrokenrationalConstitutiveLaw_TestSuite;
  }
}  // namespace CONTACT

// class implementation
class CONTACT::CONSTITUTIVELAW::BrokenrationalConstitutiveLaw_TestSuite : public CxxTest::TestSuite
{
 public:
  //! set up unit testing
  void setUp()
  {
    // initialize container for material parameters
    const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Container(1,
            INPAR::CONTACT::ConstitutiveLawType::colaw_brokenrational,
            "Brokenrational Constitutivelaw"));

    // add parameters to container
    container->Add("A", -2.);
    container->Add("B", 4.);
    container->Add("C", -0.5);
    container->Add("Offset", 0.5);

    const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> brokenrationalcoconstlaw =
        CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(container);
    coconstlaw_ = brokenrationalcoconstlaw;

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
    TS_ASSERT_DELTA(coconstlaw_->Evaluate(-2.5), -0.5, 1.e-15);

    return;
  }

  //! test member function Evaluate
  void test_EvaluateDeriv()
  {
    TS_ASSERT_DELTA(coconstlaw_->EvaluateDeriv(-2.5), 0.5, 1.e-15);
    TS_ASSERT_THROWS_ANYTHING(coconstlaw_->EvaluateDeriv(-0.25));
    return;
  }

 private:
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;
};
#endif
