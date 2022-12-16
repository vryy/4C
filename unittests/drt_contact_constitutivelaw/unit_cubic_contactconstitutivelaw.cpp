/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the power contact constitutivelaw

\level 2

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "cubic_contactconstitutivelaw.H"

namespace
{
  // class implementation
  class CubicConstitutiveLawTest : public ::testing::Test
  {
   public:
    CubicConstitutiveLawTest()
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
    }

    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;
  };


  //! test member function Evaluate
  TEST_F(CubicConstitutiveLawTest, TestEvaluate)
  {
    // gap < 0
    EXPECT_ANY_THROW(coconstlaw_->Evaluate(1.0));
    // 0< gap < offset
    EXPECT_ANY_THROW(coconstlaw_->Evaluate(-0.25));
    // offset < gap
    EXPECT_NEAR(coconstlaw_->Evaluate(-0.75), -0.8984375, 1.e-15);
  }

  //! test member function EvaluateDeriv
  TEST_F(CubicConstitutiveLawTest, TestEvaluateDeriv)
  {
    EXPECT_NEAR(coconstlaw_->EvaluateDeriv(-0.75), 4.28125, 1.e-15);
    EXPECT_ANY_THROW(coconstlaw_->EvaluateDeriv(-0.25));
  }
}  // namespace
