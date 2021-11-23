/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the linear constact constitutivelaw

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef UNIT_LINEAR_COCONSTLAW_H
#define UNIT_LINEAR_COCONSTLAW_H

#include "Unittests/src/common/unit_cxx_test_wrapper.H"

#include "../../../src/drt_contact_constitutivelaw/linear_contactconstitutivelaw.H"

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class LinearConstitutiveLaw_TestSuite;
  }
}  // namespace CONTACT

// class implementation
class CONTACT::CONSTITUTIVELAW::LinearConstitutiveLaw_TestSuite : public BACICxxTestWrapper
{
 public:
  //! set up unit testing
  void Setup() override
  {
    // initialize container for material parameters
    const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Container(
            1, INPAR::CONTACT::ConstitutiveLawType::colaw_linear, "Linear Constitutivelaw"));

    // add parameters to container
    container->Add("A", 1.5);
    container->Add("B", 0.0);
    container->Add("Offset", 0.5);

    const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> linearcoconstlaw =
        CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(container);
    coconstlaw_ = linearcoconstlaw;

    return;
  }

  //! tear down unit testing
  void TearDown() override
  {
    coconstlaw_ = Teuchos::null;
    return;
  }

  //! test member function Evaluate
  void TestEvaluate()
  {
    // gap < 0
    TS_ASSERT_THROWS_ANYTHING(coconstlaw_->Evaluate(1.0));
    // 0< gap < offset
    TS_ASSERT_THROWS_ANYTHING(coconstlaw_->Evaluate(-0.25));
    // offset < gap
    TS_ASSERT_DELTA(coconstlaw_->Evaluate(-0.75), -0.375, 1.e-15);

    return;
  }

  //! test member function Evaluate
  void TestEvaluateDeriv()
  {
    TS_ASSERT_DELTA(coconstlaw_->EvaluateDeriv(-0.75), 1.5, 1.e-15);
    TS_ASSERT_THROWS_ANYTHING(coconstlaw_->EvaluateDeriv(-0.25));
    return;
  }

 private:
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;
};
#endif
