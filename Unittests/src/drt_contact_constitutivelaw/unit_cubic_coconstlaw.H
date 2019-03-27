/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the power contact constitutivelaw

\level 2

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/
#ifndef UNIT_CUBIC_COCONSTLAW_H
#define UNIT_CUBIC_COCONSTLAW_H

#include <cxxtest/TestSuite.h>

#include "src/drt_contact_constitutivelaw/cubic_coconstlaw.H"

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class CubicConstitutiveLaw_TestSuite;
  }
}  // namespace CONTACT

// class implementation
class CONTACT::CONSTITUTIVELAW::CubicConstitutiveLaw_TestSuite : public CxxTest::TestSuite
{
 public:
  //! set up unit testing
  void setUp()
  {
    // initialize container for material parameters
    const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
        Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Container(
            1, INPAR::CONTACT::ConstitutiveLawType::colaw_cubic, "Cubic Constitutivelaw"));

    // add parameters to container
    container->Add("A", 1.5);
    container->Add("B", 2.0);
    container->Add("C", 3.0);
    container->Add("D", 0.0);
    container->Add("Offset", 0.5);

    const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> cubiccoconstlaw =
        CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(container);
    coconstlaw_ = cubiccoconstlaw;

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
    TS_ASSERT_DELTA(coconstlaw_->Evaluate(-0.75), -0.8984375, 1.e-15);

    return;
  }

  //! test member function Evaluate
  void test_EvaluateDeriv()
  {
    TS_ASSERT_DELTA(coconstlaw_->EvaluateDeriv(-0.75), 4.28125, 1.e-15);
    return;
  }

 private:
  Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;
};
#endif
